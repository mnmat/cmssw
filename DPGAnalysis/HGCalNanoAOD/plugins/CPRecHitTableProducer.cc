#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "PhysicsTools/NanoAOD/interface/SimpleFlatTableProducer.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/TICLGeomTools.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Validation/HGCalValidation/interface/CaloParticleSelector.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticleFwd.h"
#include "SimDataFormats/CaloAnalysis/interface/SimClusterFwd.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"

class CPRecHitTableProducer : public edm::stream::EDProducer<edm::stream::WatchRuns> {
public:
  CPRecHitTableProducer(edm::ParameterSet const& params)
      : skipNonExistingSrc_(params.getParameter<bool>("skipNonExistingSrc")),
        hitTableName_(params.getParameter<std::string>("hitTableName")),
        trackTableName_(params.getParameter<std::string>("trackTableName")),
        precision_(params.getParameter<int>("precision")),
        caloParticlesToken_(consumes<std::vector<CaloParticle> >(params.getParameter<edm::InputTag>("caloParticles"))), 
        hgcalRecHitsEEToken_(consumes<HGCRecHitCollection>(params.getParameter<edm::InputTag>("hgcalRecHitsEE"))),
        hgcalRecHitsFHToken_(consumes<HGCRecHitCollection>(params.getParameter<edm::InputTag>("hgcalRecHitsFH"))),
        hgcalRecHitsBHToken_(consumes<HGCRecHitCollection>(params.getParameter<edm::InputTag>("hgcalRecHitsBH"))),
        simVerticesToken_(consumes<std::vector<SimVertex>>(params.getParameter<edm::InputTag>("simVertices"))),
        ticlGeomToken_(esConsumes<TICLGeomHost, CaloGeometryRecord, edm::Transition::BeginRun>(edm::ESInputTag("", ""))),
        ticlGeomLookupToken_(
            esConsumes<TICLGeomLookupHost, CaloGeometryRecord, edm::Transition::BeginRun>(edm::ESInputTag("", ""))),
        ticlGeomLayersToken_(
            esConsumes<TICLGeomLayersHost, CaloGeometryRecord, edm::Transition::BeginRun>(edm::ESInputTag("", ""))) {
    produces<nanoaod::FlatTable>(hitTableName_);
    produces<nanoaod::FlatTable>(trackTableName_);
  }

  void fillHitMap(std::map<DetId, HGCRecHit>& hitMap,
    const HGCRecHitCollection& rechitsEE,
    const HGCRecHitCollection& rechitsFH,
    const HGCRecHitCollection& rechitsBH) const {
    hitMap.clear();
    int counter = 0;
    for (const auto& hit : rechitsEE) {
      hitMap.emplace(hit.detid(), hit);
      counter++;
    }
  
    for (const auto& hit : rechitsFH) {
      hitMap.emplace(hit.detid(), hit);
      counter++;
    }
  
    for (const auto& hit : rechitsBH) {
      hitMap.emplace(hit.detid(), hit);
      counter++;
    }
  } // end of EfficiencyStudies::fillHitMap


  void beginRun(edm::Run const&, edm::EventSetup const& es) override {
    rhtools_.setGeometry(es.getData(ticlGeomToken_), es.getData(ticlGeomLookupToken_), es.getData(ticlGeomLayersToken_));
  }

  void produce(edm::Event& iEvent, const edm::EventSetup& es) override {
    //Layer Clusters time value map

    edm::Handle<std::vector<SimVertex>> simVerticesHandle;
    iEvent.getByToken(simVerticesToken_, simVerticesHandle);
    std::vector<SimVertex> const& simVertices = *simVerticesHandle;

    edm::Handle<HGCRecHitCollection> recHitHandleEE;  
    iEvent.getByToken(hgcalRecHitsEEToken_, recHitHandleEE);
  
    edm::Handle<HGCRecHitCollection> recHitHandleFH;
    iEvent.getByToken(hgcalRecHitsFHToken_, recHitHandleFH);
  
    edm::Handle<HGCRecHitCollection> recHitHandleBH;
    iEvent.getByToken(hgcalRecHitsBHToken_, recHitHandleBH);

    std::map<DetId, HGCRecHit> hitMap;
    fillHitMap(hitMap, *recHitHandleEE, *recHitHandleFH, *recHitHandleBH);


    auto cps_h = iEvent.getHandle(caloParticlesToken_);

    const auto ncps = cps_h.isValid() ? cps_h->size() : 0;


    // Define CP selector
    CaloParticleSelector cpSelector_ = CaloParticleSelector(0, 100, -3, 3,120,280,0,1000000,true,false,false,false,false,{13},-3.2,3.2);

    static constexpr float default_value = std::numeric_limits<float>::quiet_NaN();

    std::vector<int> detids;
    std::vector<float> xs;
    std::vector<float> ys;
    std::vector<float> zs;
    std::vector<uint16_t> layers;
    std::vector<uint16_t> nHits(ncps, 0);
    
    if (cps_h.isValid() || !(this->skipNonExistingSrc_)) {
      const auto& cps = *cps_h;
      for (size_t cpIndex = 0; cpIndex < ncps; ++cpIndex) {        // do something with track parameters, e.g, plot the charge.
        const CaloParticle& cp = cps[cpIndex]; 
    
        //Select only Signal
        if (!cpSelector_(cp,simVertices)) continue;
        auto pid = cp.particleId();
        const SimClusterRefVector& simclusters = cp.simClusters();
    
        int simClusterId = 0;
        for (const auto& it_simc : simclusters){
          const SimCluster& simc = (*(it_simc));
          const auto& sc_haf = simc.hits_and_fractions();
    
          for (const auto& it_sc_haf : sc_haf){
    
            DetId detid_ = (it_sc_haf.first);

            std::map<DetId, HGCRecHit>::const_iterator itcheck = hitMap.find(detid_);

            if (itcheck != hitMap.end()){
              auto gp = rhtools_.getPosition(detid_);

              xs.push_back(gp.x());
              ys.push_back(gp.y());
              zs.push_back(gp.z());
              layers.push_back(rhtools_.getLayerWithOffset(detid_));
              detids.push_back(detid_);
              ++nHits[cpIndex];
              // process the recHit as needed
            }
          }
        }
      } 
    }

    auto counts = std::make_unique<nanoaod::FlatTable>(ncps, trackTableName_, /*singleton*/ false, /*extension*/ false);
    counts->addColumn<uint16_t>("nHits", nHits, "number of valid HGCAL rechits on the CP (see " + hitTableName_ + ")");

    auto cpRecHitTable =
        std::make_unique<nanoaod::FlatTable>(xs.size(), hitTableName_, /*singleton*/ false, /*extension*/ false);
    cpRecHitTable->addColumn<int>("detid", detids, "CP RecHit detid", precision_);
    cpRecHitTable->addColumn<float>("x", xs, "HGC Track RecHit x position [cm]", precision_);
    cpRecHitTable->addColumn<float>("y", ys, "HGC Track RecHit y position [cm]", precision_);
    cpRecHitTable->addColumn<float>("z", zs, "HGC Track RecHit z position [cm]", precision_);
    cpRecHitTable->addColumn<uint16_t>("layer", layers, "HGC Track RecHit layer number", precision_);
    iEvent.put(std::move(cpRecHitTable), hitTableName_);
    iEvent.put(std::move(counts), trackTableName_);
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<bool>("skipNonExistingSrc", false)
        ->setComment("whether or not to skip producing the table on absent input product");
    desc.add<std::string>("hitTableName", "HGCTrackRecHits")->setComment("name of the flat table ouput");
    desc.add<std::string>("trackTableName", "HGCTrackRecHits_links")->setComment("name of the flat table ouput");
    desc.add<edm::InputTag>("caloParticles", edm::InputTag("mix", "MergedCaloTruth"))->setComment("input collection of CaloParticles");
    desc.add<edm::InputTag>("hgcalRecHitsEE", edm::InputTag("HGCalRecHit", "HGCEERecHits"));
    desc.add<edm::InputTag>("hgcalRecHitsFH", edm::InputTag("HGCalRecHit", "HGCHEFRecHits"));
    desc.add<edm::InputTag>("hgcalRecHitsBH", edm::InputTag("HGCalRecHit", "HGCHEBRecHits"));
    desc.add<edm::InputTag>("simVertices", edm::InputTag("mix", "MergedTruth"))->setComment("input collection of SimVertices");
    desc.add<int>("precision", 7);
    descriptions.addWithDefaultLabel(desc);
  }

private:
  const bool skipNonExistingSrc_;
  const std::string hitTableName_;
  const std::string trackTableName_;
  const unsigned int precision_;
  edm::EDGetTokenT<std::vector<CaloParticle>> caloParticlesToken_;
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsEEToken_; 
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsFHToken_;
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsBHToken_;
  edm::EDGetTokenT<std::vector<SimVertex>> simVerticesToken_;

  edm::ESGetToken<TICLGeomHost, CaloGeometryRecord> ticlGeomToken_;
  edm::ESGetToken<TICLGeomLookupHost, CaloGeometryRecord> ticlGeomLookupToken_;
  edm::ESGetToken<TICLGeomLayersHost, CaloGeometryRecord> ticlGeomLayersToken_;
  ticlgeom::Tools rhtools_;

};

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(CPRecHitTableProducer);