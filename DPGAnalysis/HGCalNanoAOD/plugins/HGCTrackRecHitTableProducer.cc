#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "PhysicsTools/NanoAOD/interface/SimpleFlatTableProducer.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/TICLGeomTools.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"

class HGCTrackRecHitTableProducer : public edm::stream::EDProducer<edm::stream::WatchRuns> {
public:
  HGCTrackRecHitTableProducer(edm::ParameterSet const& params)
      : skipNonExistingSrc_(params.getParameter<bool>("skipNonExistingSrc")),
        hitTableName_(params.getParameter<std::string>("hitTableName")),
        trackTableName_(params.getParameter<std::string>("trackTableName")),
        precision_(params.getParameter<int>("precision")),
        hgcTracksToken_(consumes<edm::View<reco::Track>>(params.getParameter<edm::InputTag>("hgcTracks"))),
        hgcalRecHitsEEToken_(consumes<HGCRecHitCollection>(params.getParameter<edm::InputTag>("hgcalRecHitsEE"))),
        hgcalRecHitsFHToken_(consumes<HGCRecHitCollection>(params.getParameter<edm::InputTag>("hgcalRecHitsFH"))),
        hgcalRecHitsBHToken_(consumes<HGCRecHitCollection>(params.getParameter<edm::InputTag>("hgcalRecHitsBH"))),
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

    edm::Handle<HGCRecHitCollection> recHitHandleEE;  
    iEvent.getByToken(hgcalRecHitsEEToken_, recHitHandleEE);
  
    edm::Handle<HGCRecHitCollection> recHitHandleFH;
    iEvent.getByToken(hgcalRecHitsFHToken_, recHitHandleFH);
  
    edm::Handle<HGCRecHitCollection> recHitHandleBH;
    iEvent.getByToken(hgcalRecHitsBHToken_, recHitHandleBH);

    std::map<DetId, HGCRecHit> hitMap;
    fillHitMap(hitMap, *recHitHandleEE, *recHitHandleFH, *recHitHandleBH);

    auto hgcTracks_h = iEvent.getHandle(hgcTracksToken_);
    const auto nhgcTracks = hgcTracks_h.isValid() ? hgcTracks_h->size() : 0;

    static constexpr float default_value = std::numeric_limits<float>::quiet_NaN();

    std::vector<int> detids;
    std::vector<float> xs;
    std::vector<float> ys;
    std::vector<float> zs;
    std::vector<uint16_t> layers;
    std::vector<uint16_t> nHits(nhgcTracks, 0);
    
    if (hgcTracks_h.isValid() || !(this->skipNonExistingSrc_)) {
      const auto& tracks = *hgcTracks_h;
      for (size_t tkIndex = 0; tkIndex < nhgcTracks; ++tkIndex) {
        const auto& track = tracks[tkIndex];
        for (auto it = track.recHitsBegin(); it != track.recHitsEnd(); ++it) {
          const auto& recHit = **it;
          auto detid = recHit.rawId();
          std::map<DetId, HGCRecHit>::const_iterator itcheck = hitMap.find(detid);

          if (itcheck != hitMap.end()){
            auto gp = rhtools_.getPosition(detid);

            xs.push_back(gp.x());
            ys.push_back(gp.y());
            zs.push_back(gp.z());
            layers.push_back(rhtools_.getLayerWithOffset(detid));
            detids.push_back(detid);
            ++nHits[tkIndex];
          }
          // process the recHit as needed
        }
      }
    } 

    auto counts = std::make_unique<nanoaod::FlatTable>(nhgcTracks, trackTableName_, /*singleton*/ false, /*extension*/ false);
    counts->addColumn<uint16_t>("nHits", nHits, "number of valid HGCAL rechits on the track (see " + hitTableName_ + ")");

    auto hgcTrackRecHitTable =
        std::make_unique<nanoaod::FlatTable>(xs.size(), hitTableName_, /*singleton*/ false, /*extension*/ false);
    hgcTrackRecHitTable->addColumn<int>("detid", detids, "HGC Track RecHit detid", precision_);
    hgcTrackRecHitTable->addColumn<float>("x", xs, "HGC Track RecHit x position [cm]", precision_);
    hgcTrackRecHitTable->addColumn<float>("y", ys, "HGC Track RecHit y position [cm]", precision_);
    hgcTrackRecHitTable->addColumn<float>("z", zs, "HGC Track RecHit z position [cm]", precision_);
    hgcTrackRecHitTable->addColumn<uint16_t>("layer", layers, "HGC Track RecHit layer number", precision_);
    iEvent.put(std::move(hgcTrackRecHitTable), hitTableName_);
    iEvent.put(std::move(counts), trackTableName_);
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<bool>("skipNonExistingSrc", false)
        ->setComment("whether or not to skip producing the table on absent input product");
    desc.add<std::string>("hitTableName", "HGCTrackRecHits")->setComment("name of the flat table ouput");
    desc.add<std::string>("trackTableName", "HGCTrackRecHits_links")->setComment("name of the flat table ouput");
    desc.add<edm::InputTag>("hgcTracks", edm::InputTag("ticlTrackstersKalmanFilter","HGCTracks","RECO"))
        ->setComment("input collection of HGC tracks");
    desc.add<edm::InputTag>("hgcalRecHitsEE", edm::InputTag("HGCalRecHit", "HGCEERecHits"));
    desc.add<edm::InputTag>("hgcalRecHitsFH", edm::InputTag("HGCalRecHit", "HGCHEFRecHits"));
    desc.add<edm::InputTag>("hgcalRecHitsBH", edm::InputTag("HGCalRecHit", "HGCHEBRecHits"));
    desc.add<int>("precision", 7);
    descriptions.addWithDefaultLabel(desc);
  }

private:
  const bool skipNonExistingSrc_;
  const std::string hitTableName_;
  const std::string trackTableName_;
  const unsigned int precision_;
  edm::EDGetTokenT<edm::View<reco::Track>> hgcTracksToken_;
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsEEToken_; 
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsFHToken_;
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsBHToken_;

  edm::ESGetToken<TICLGeomHost, CaloGeometryRecord> ticlGeomToken_;
  edm::ESGetToken<TICLGeomLookupHost, CaloGeometryRecord> ticlGeomLookupToken_;
  edm::ESGetToken<TICLGeomLayersHost, CaloGeometryRecord> ticlGeomLayersToken_;
  ticlgeom::Tools rhtools_;

};

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(HGCTrackRecHitTableProducer);