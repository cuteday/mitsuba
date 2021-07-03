/* BDPM modified on PM */

#include <mitsuba/core/plugin.h>
#include <mitsuba/render/common.h>
#include <mitsuba/render/gatherproc.h>
#include "bre.h"

MTS_NAMESPACE_BEGIN

class BDPMIntegrator : public SamplingIntegrator {
public:
    BDPMIntegrator(const Properties &props) : SamplingIntegrator(props),
          m_parentIntegrator(NULL) {
        /* Depth to start using russian roulette when tracing photons */
        m_rrDepth = props.getInteger("rrDepth", 5);
        /* Longest visualized path length (\c -1 = infinite).
           A value of \c 1 will visualize only directly visible light sources.
           \c 2 will lead to single-bounce (direct-only) illumination, and so on. */
        m_maxDepth = props.getInteger("maxDepth", -1);
        /* Longest depth of camera subpath. */
        m_maxCameraDepth = props.getInteger("maxCameraDepth", -1);
        /* Maxium diffuse bounces in camera subpath. */
        m_maxCameraDiffuseBounce = props.getInteger("maxCameraDiffuseBounce", -1);
        /* Maxium glossy bounces in camera subpath. */
        m_maxCameraGlossyBounce = props.getInteger("maxCameraGlossyBounce", -1);
        /* Granularity of photon tracing work units (in shot particles, 0 => decide automatically) */
        m_granularity = props.getInteger("granularity", 0);
        /* Number of photons to collect for the global photon map */
        m_globalPhotons = props.getSize("globalPhotons", 250000);
        /* Max. radius of lookups in the global photon map (relative to the scene size) */
        m_globalLookupRadiusRel = props.getFloat("globalLookupRadius", 0.05f);
        /* Minimum amount of photons to consider a photon map lookup valid */
        int lookupSize = props.getInteger("lookupSize", 120);
        /* Minimum amount of photons to consider a volumetric photon map lookup valid */
        m_globalLookupSize = props.getInteger("globalLookupSize", lookupSize);
        /* Should photon gathering steps exclusively run on the local machine? */
        m_gatherLocally = props.getBoolean("gatherLocally", true);
        /* Indicates if the gathering steps should be canceled if not enough photons are generated. */
        m_autoCancelGathering = props.getBoolean("autoCancelGathering", true);
        /* When this flag is set to true, contributions from directly
         * visible emitters will not be included in the rendered image */
        m_hideEmitters = props.getBoolean("hideEmitters", false);

        if (m_maxDepth == 0 || m_maxCameraDepth == 0) {
            Log(EError, "maxDepth must be greater than zero!");
        } else {
            /**
             * An infinite depth is currently not supported, since
             * the photon tracing step uses a Halton sequence
             * that is based on a finite-sized prime number table
             */
            if (m_maxDepth == -1) 
                m_maxDepth = 128;
            if (m_maxCameraDepth == -1)
                m_maxCameraDepth = 128;
            if (m_maxCameraDiffuseBounce == -1)
                m_maxCameraDiffuseBounce = 128;
            if (m_maxCameraGlossyBounce == -1)
                m_maxCameraGlossyBounce = 128;
            Log(EInfo, "maxDepth set to %d", m_maxDepth);
            Log(EInfo, "maxCameraDepth set to %d", m_maxCameraDepth);
            Log(EInfo, "maxCameraDiffuseBounce set to %d", m_maxCameraDiffuseBounce);
            Log(EInfo, "maxCameraGlossyBounce set to %d", m_maxCameraGlossyBounce);
        }

        m_globalPhotonMapID = m_breID = 0;
    }

    /// Unserialize from a binary data stream
    BDPMIntegrator(Stream *stream, InstanceManager *manager)
     : SamplingIntegrator(stream, manager), m_parentIntegrator(NULL) {
        m_maxDepth = stream->readInt();
        m_rrDepth = stream->readInt();
        m_globalPhotons = stream->readSize();
        m_globalLookupRadius = stream->readFloat();
        m_globalLookupSize = stream->readInt();
        m_gatherLocally = stream->readBool();
        m_autoCancelGathering = stream->readBool();
        m_hideEmitters = stream->readBool();
        m_globalPhotonMapID = m_breID = 0;
        configure();
    }

    virtual ~BDPMIntegrator() {
        ref<Scheduler> sched = Scheduler::getInstance();
        if (m_globalPhotonMapID)
            sched->unregisterResource(m_globalPhotonMapID);
        if (m_breID)
            sched->unregisterResource(m_breID);
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        SamplingIntegrator::serialize(stream, manager);
        stream->writeInt(m_maxDepth);
        stream->writeInt(m_rrDepth);
        stream->writeSize(m_globalPhotons);
        stream->writeFloat(m_globalLookupRadius);
        stream->writeInt(m_globalLookupSize);
        stream->writeBool(m_gatherLocally);
        stream->writeBool(m_autoCancelGathering);
        stream->writeBool(m_hideEmitters);
    }

    /// Configure the sampler for a specified amount of direct illumination samples
    void configureSampler(const Scene *scene, Sampler *sampler) {
        SamplingIntegrator::configureSampler(scene, sampler);
    }

    void configure() {}

    bool preprocess(const Scene *scene, RenderQueue *queue, const RenderJob *job,
            int sceneResID, int sensorResID, int samplerResID) {
        SamplingIntegrator::preprocess(scene, queue, job, sceneResID, sensorResID, samplerResID);
        /* Create a deterministic sampler for the photon gathering step */
        ref<Scheduler> sched = Scheduler::getInstance();
        ref<Sampler> sampler = static_cast<Sampler *> (PluginManager::getInstance()->
            createObject(MTS_CLASS(Sampler), Properties("halton")));
        /* Create a sampler instance for every core */
        std::vector<SerializableObject *> samplers(sched->getCoreCount());
        for (size_t i=0; i<sched->getCoreCount(); ++i) {
            ref<Sampler> clonedSampler = sampler->clone();
            clonedSampler->incRef();
            samplers[i] = clonedSampler.get();
        }
        int qmcSamplerID = sched->registerMultiResource(samplers);
        for (size_t i=0; i<samplers.size(); ++i)
            samplers[i]->decRef();

        const ref_vector<Medium> &media = scene->getMedia();

        for (ref_vector<Medium>::const_iterator it = media.begin(); it != media.end(); ++it) {
            Log(EWarn, "An media found!");
            if (!(*it)->isHomogeneous())
                Log(EError, "Inhomogeneous media are currently not supported by the photon mapper!");
        }

        if (m_globalPhotonMap.get() == NULL && m_globalPhotons > 0) {
            /* Generate the global photon map */
            ref<GatherPhotonProcess> proc = new GatherPhotonProcess(
                GatherPhotonProcess::EAllSurfacePhotons, m_globalPhotons,
                m_granularity, m_maxDepth-1, m_rrDepth, m_gatherLocally,
                m_autoCancelGathering, job);

            proc->bindResource("scene", sceneResID);
            proc->bindResource("sensor", sensorResID);
            proc->bindResource("sampler", qmcSamplerID);

            m_proc = proc;
            sched->schedule(proc);
            sched->wait(proc);
            m_proc = NULL;

            if (proc->getReturnStatus() != ParallelProcess::ESuccess)
                return false;

            ref<PhotonMap> globalPhotonMap = proc->getPhotonMap();
            if (globalPhotonMap->isFull()) {
                Log(EDebug, "Global photon map full. Shot " SIZE_T_FMT " particles, excess photons due to parallelism: "
                    SIZE_T_FMT, proc->getShotParticles(), proc->getExcessPhotons());

                m_globalPhotonMap = globalPhotonMap;
                m_globalPhotonMap->setScaleFactor(1 / (Float) proc->getShotParticles());
                m_globalPhotonMap->build();
                m_globalPhotonMapID = sched->registerResource(m_globalPhotonMap);
            }
        }

        
        /* Adapt to scene extents */
        m_globalLookupRadius = m_globalLookupRadiusRel * scene->getBSphere().radius;
        sched->unregisterResource(qmcSamplerID);
        return true;
    }

    void setParent(ConfigurableObject *parent) {
        if (parent->getClass()->derivesFrom(MTS_CLASS(SamplingIntegrator)))
            m_parentIntegrator = static_cast<SamplingIntegrator *>(parent);
        else
            m_parentIntegrator = this;
    }

    /// Specify globally shared resources
    void bindUsedResources(ParallelProcess *proc) const {
        if (m_globalPhotonMap.get())
            proc->bindResource("globalPhotonMap", m_globalPhotonMapID);
        if (m_bre.get())
            proc->bindResource("bre", m_breID);
    }

    /// Connect to globally shared resources
    void wakeup(ConfigurableObject *parent, std::map<std::string, SerializableObject *> &params) {
        if (!m_globalPhotonMap.get() && params.find("globalPhotonMap") != params.end())
            m_globalPhotonMap = static_cast<PhotonMap *>(params["globalPhotonMap"]);
        if (!m_bre.get() && params.find("bre") != params.end())
            m_bre = static_cast<BeamRadianceEstimator *>(params["bre"]);
        if (parent && parent->getClass()->derivesFrom(MTS_CLASS(SamplingIntegrator)))
            m_parentIntegrator = static_cast<SamplingIntegrator *>(parent);
        else
            m_parentIntegrator = this;
    }

    void cancel() {
        SamplingIntegrator::cancel();
        if (m_proc)
            Scheduler::getInstance()->cancel(m_proc);
    }

    Spectrum Li(const RayDifferential &ray, RadianceQueryRecord &rRec) const {
        Log(EDebug, "Entering Li depth %d", rRec.depth);

        Spectrum LiSurf(0.0f);
        Intersection &its = rRec.its;
        const Scene *scene = rRec.scene;

		Spectrum throughput(1.0f);

        /* Perform the first ray intersection (or ignore if the intersection has already been provided). */
        rRec.rayIntersect(ray);

		while ((rRec.depth <= m_maxDepth && m_maxDepth > 0 ) || (rRec.depth <= m_maxCameraDepth && m_maxCameraDepth > 0 )){
			if (!its.isValid()) {
				/* If no intersection could be found, possibly return
				attenuated radiance from a background luminaire */
				if ((rRec.type & RadianceQueryRecord::EEmittedRadiance) && !m_hideEmitters)
					LiSurf += throughput *scene->evalEnvironment(ray);
				break;
			}

			/* Possibly include emitted radiance if requested */
			if (its.isEmitter() && (rRec.type & RadianceQueryRecord::EEmittedRadiance) && !m_hideEmitters)
				LiSurf += throughput *its.Le(-ray.d);

			/* Subsurface */
			if (its.hasSubsurface() && (rRec.type & RadianceQueryRecord::ESubsurfaceRadiance))
				LiSurf += throughput *its.LoSub(scene, rRec.sampler, -ray.d, rRec.depth);

			const BSDF *bsdf = its.getBSDF(ray);
			unsigned int bsdfType = bsdf->getType() & BSDF::EAll;

			/* Irradiance cache query -> treat as diffuse */
			if (bsdfType & BSDF::ESmooth) {
				/* Estimate radiance using photon map on diffuse surfaces */
				int maxDepth = m_maxDepth == -1 ? INT_MAX : (m_maxDepth-rRec.depth);
				if (rRec.type & RadianceQueryRecord::EIndirectSurfaceRadiance && m_globalPhotonMap.get()){
					LiSurf += throughput * m_globalPhotonMap->estimateRadianceBDPM(its, m_globalLookupRadius, m_globalLookupSize,
						maxDepth, rRec.pathProb, rRec.invPdf, m_rrDepth, 0.8f);
				}
			}

			/* Determine if the ray is absorbed by the RR strategy. */
			Float judgeRR = rRec.nextSample1D();
			Float probRR = rRec.depth > m_rrDepth ? 0.8f : 1.0f;
			if(judgeRR > probRR){
				break;
			}

			/*               BSDF sampling               */
			Point2 sample = rRec.nextSample2D();

			/* Sample BSDF * cos(theta) */
			BSDFSamplingRecord bRec(its, rRec.sampler, ERadiance);

			Spectrum bsdfVal = bsdf->sample(bRec, sample);
			//BSDFSamplingRecord bInvRec(its, bRec.wo, bRec.wi);
			BSDFSamplingRecord bInvRec = bRec;
			bInvRec.wi = bRec.wo, bInvRec.wo = bRec.wi;

			Float bsdfPdf, invBsdfPdf;
			bsdfPdf = bsdf->pdf(bRec, bRec.sampledType & BSDF::EDelta ? EDiscrete : ESolidAngle);
			invBsdfPdf = bsdf->pdf(bInvRec, bRec.sampledType & BSDF::EDelta ? EDiscrete : ESolidAngle);

			if (bsdfVal.isZero())
				break;

			/* Trace a ray in this direction */
			RayDifferential bsdfRay(its.p, its.toWorld(bRec.wo), ray.time);
			scene->rayIntersect(bsdfRay, its);

			if (bRec.sampledType & BSDF::EDelta){
				if (++rRec.glossyBounce >= m_maxCameraGlossyBounce)
					break;
			}
			else if (bRec.sampledType & BSDF::ESmooth){
				if (++rRec.diffuseBounce >= m_maxCameraDiffuseBounce)
					break;
			}

			throughput *= bsdfVal / probRR;

			// [for bdpm] recursively added path prob records
			rRec.depth++;
			rRec.pathProb.push_back(rRec.pathProb.back() * bsdfPdf);
			rRec.invPdf.push_back(invBsdfPdf);

			//Log(EInfo, "Current depth %d, pathProb length %d, invPdf length %d", rRec2.depth, rRec2.pathProb.size(), rRec2.invPdf.size());
		
		}
        //Log(EInfo, "Current LiSurf: %f", LiSurf);

        return LiSurf;
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "BDPMIntegrator[" << endl
            << "  maxDepth = " << m_maxDepth << "," << endl
            << "  rrDepth = " << m_rrDepth << "," << endl
            << "  globalPhotons = " << m_globalPhotons << "," << endl
            << "  gatherLocally = " << m_gatherLocally << "," << endl
            << "  globalLookupRadius = " << m_globalLookupRadius << "," << endl
            << "  globalLookupSize = " << m_globalLookupSize << "," << endl
            << "]";
        return oss.str();
    }

    /* Power(2) heuristic for multi-importance sampling. */
    inline Float miWeight(Float pdfA, Float pdfB) const {
        pdfA *= pdfA; pdfB *= pdfB;
        return pdfA / (pdfA + pdfB);
    }

    MTS_DECLARE_CLASS()
private:
    ref<PhotonMap> m_globalPhotonMap;
    ref<BeamRadianceEstimator> m_bre;
    ref<ParallelProcess> m_proc;
    SamplingIntegrator *m_parentIntegrator;
    int m_globalPhotonMapID, m_breID;
    size_t m_globalPhotons;
    int m_globalLookupSize;
    Float m_globalLookupRadiusRel, m_globalLookupRadius;
    int m_granularity;
    int m_rrDepth, m_maxDepth;
    int m_maxCameraDepth, m_maxCameraDiffuseBounce, m_maxCameraGlossyBounce;
    bool m_gatherLocally, m_autoCancelGathering;
    bool m_hideEmitters;
};

MTS_IMPLEMENT_CLASS_S(BDPMIntegrator, false, SamplingIntegrator)
MTS_EXPORT_PLUGIN(BDPMIntegrator, "Photon map integrator");
MTS_NAMESPACE_END
