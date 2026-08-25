// Author: Mark Matthewman - mark.matthewman@cern.ch
// Date: 07/2026

#include <algorithm>
#include <memory>
#include <unordered_map>
#include <vector>

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/PluginDescription.h"

#include "DataFormats/HGCalReco/interface/KFHit.h"
#include "DataFormats/HGCalReco/interface/TICLLayerTile.h"
#include "DataFormats/HGCalReco/interface/TICLSeedingRegion.h"

#include "RecoHGCal/TICL/plugins/HGCTrackingPluginFactory.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/PatternTools/interface/TempTrajectory.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"


class HGCTracksProducer
    : public edm::stream::EDProducer<edm::stream::WatchRuns>{
public:
  explicit HGCTracksProducer(const edm::ParameterSet&);
  ~HGCTracksProducer() override {}

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  void produce(edm::Event&, const edm::EventSetup&) override;

  void beginRun(const edm::Run&, const edm::EventSetup& es) override {
    const auto& geom = es.getData(geometry_token_);
    rhtools_.setGeometry(geom);
  
    myAlgo_->setGeometry(rhtools_);
  }
  
private:
  std::unique_ptr<ticl::HGCTrackingAlgoBaseT<TICLLayerTiles>> myAlgo_;

  edm::EDGetTokenT<TICLLayerTiles> rechit_tiles_token_;

  const edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geometry_token_;
  hgcal::RecHitTools rhtools_;
  const std::string itername_;
  const edm::EDGetTokenT<std::vector<TICLSeedingRegion>> seeding_regions_token_;

};

DEFINE_FWK_MODULE(HGCTracksProducer);

HGCTracksProducer::HGCTracksProducer(const edm::ParameterSet& ps)
    : geometry_token_(esConsumes<CaloGeometry, CaloGeometryRecord, edm::Transition::BeginRun>()),
      itername_(ps.getParameter<std::string>("itername")),
      seeding_regions_token_(
        consumes<std::vector<TICLSeedingRegion>>(ps.getParameter<edm::InputTag>("seeding_regions"))
      ){
  const auto plugin = ps.getParameter<std::string>("hgcTrackingBy");
  const auto pluginPSet = ps.getParameter<edm::ParameterSet>("pluginHGCTrackingBy" + plugin);

  myAlgo_ = HGCTrackingFactory::get()->create(plugin, pluginPSet, consumesCollector());
  rechit_tiles_token_ = consumes<TICLLayerTiles>(ps.getParameter<edm::InputTag>("rechit_tiles"));


  produces<std::vector<KFHit>>("KFHits");
  produces<std::vector<reco::Track>>("HGCALTracks");
  produces<std::vector<reco::TrackExtra>>("HGCALTrackExtras");
  produces<TrackingRecHitCollection>("HGCALTrackingRecHitCollection");
}

void HGCTracksProducer::produce(edm::Event& evt, const edm::EventSetup& es) {
  auto kfhits = std::make_unique<std::vector<KFHit>>();
  auto tracks = std::make_unique<std::vector<reco::Track>>();
  auto trackExtras = std::make_unique<std::vector<reco::TrackExtra>>();
  auto trackingRecHitCollection = std::make_unique<TrackingRecHitCollection>();

  const auto& tiles = evt.get(rechit_tiles_token_);
  const auto& seeding_regions = evt.get(seeding_regions_token_);

  if (!seeding_regions.empty()){
    const typename ticl::HGCTrackingAlgoBaseT<TICLLayerTiles>::Inputs input(evt, es, tiles, seeding_regions);
    myAlgo_->makeTrajectories(input, *kfhits, *tracks, *trackExtras, *trackingRecHitCollection);
  }

  evt.put(std::move(kfhits), "KFHits");
  evt.put(std::move(tracks), "HGCALTracks");
  evt.put(std::move(trackExtras), "HGCALTrackExtras");
  evt.put(std::move(trackingRecHitCollection), "HGCALTrackingRecHitCollection");
}

void HGCTracksProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::string>("detector", "HGCAL");
  desc.add<edm::InputTag>("rechit_tiles", edm::InputTag("hgcalLayerClustersTiles"));
  desc.add<edm::InputTag>("seeding_regions", edm::InputTag("ticlSeedingRegionProducer"));
  desc.add<std::string>("itername", "HGCTracks");
  desc.add<std::string>("hgcTrackingBy", "KalmanFilter");

  // HGCTracking plugins
  edm::ParameterSetDescription pluginDescKF;
  pluginDescKF.addNode(edm::PluginDescription<HGCTrackingFactory>("type", "KalmanFilter", true));
  desc.add<edm::ParameterSetDescription>("pluginHGCTrackingByKalmanFilter", pluginDescKF);

  descriptions.add("hgcTracksProducer", desc);
}