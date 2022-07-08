// Author: Marco Rovere, marco.rovere@cern.ch
// Date: 05/2019
//
#include <memory>  // unique_ptr

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/ESGetToken.h"

#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/HGCalReco/interface/TICLLayerTile.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/TICLGeomTools.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"


struct TileParameters{
  int layer;
  double eta, phi;
};

class TICLLayerTileProducer : public edm::stream::EDProducer<edm::stream::WatchRuns> {
  public:
  explicit TICLLayerTileProducer(const edm::ParameterSet &ps);
  ~TICLLayerTileProducer() override {}
  void beginRun(edm::Run const &, edm::EventSetup const &) override;
  void produce(edm::Event &, const edm::EventSetup &) override;
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  void produceLCTiles(edm::Event &evt, TICLLayerTiles* result, TICLLayerTilesBarrel* resultBarrel, TICLLayerTilesHFNose* resultHFNose);
  void produceRecHitTiles(edm::Event &evt, TICLLayerTiles* result);
  TileParameters getTileParameters(const HGCRecHit& hit); 

  template<typename T, typename U> void fillTiles(T& results, 
      const U& objects);

  edm::EDGetTokenT<std::vector<reco::CaloCluster>> clusters_token_;
  edm::EDGetTokenT<std::vector<reco::CaloCluster>> clusters_HFNose_token_;
  edm::ESGetToken<TICLGeomHost, CaloGeometryRecord> ticlGeomToken_;
  edm::ESGetToken<TICLGeomLookupHost, CaloGeometryRecord> ticlGeomLookupToken_;
  edm::ESGetToken<TICLGeomLayersHost, CaloGeometryRecord> ticlGeomLayersToken_;
  ticlgeom::Tools rhtools_;
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsEEToken_; 
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsFHToken_;
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsBHToken_;  
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsHFNoseToken_;

  std::string detector_;
  bool isLC_;
  bool doNose_;
  bool doBarrel_;

  int offset_ = 0;
};

void TICLLayerTileProducer::beginRun(edm::Run const &, edm::EventSetup const &es) {
  rhtools_.setGeometry(es.getData(ticlGeomToken_), es.getData(ticlGeomLookupToken_), es.getData(ticlGeomLayersToken_));
}

void TICLLayerTileProducer::produceLCTiles(edm::Event &evt,
  TICLLayerTiles* result,
  TICLLayerTilesBarrel* resultBarrel,
  TICLLayerTilesHFNose* resultHFNose) {
  edm::Handle<std::vector<reco::CaloCluster>> cluster_h;
  doNose_ ? evt.getByToken(clusters_HFNose_token_, cluster_h)
  : evt.getByToken(clusters_token_, cluster_h);

  int lcId = offset_;
  for (auto const &lc : *cluster_h) {
    const auto firstHitDetId = lc.hitsAndFractions()[0].first;
    int layer = rhtools_.getLayerWithOffset(firstHitDetId);
    bool isBarrelLC = rhtools_.isBarrel(firstHitDetId);
    if (!isBarrelLC) {
      layer += rhtools_.lastLayer(doNose_) * ((rhtools_.zside(firstHitDetId) + 1) >> 1) - 1;
    }
    assert(layer >= 0);

    if (doNose_) {
      resultHFNose->fill(layer, lc.eta(), lc.phi(), lcId);
      LogDebug("TICLLayerTileProducer") << "Adding layerClusterId: " << lcId << " into bin [eta,phi]: [ "
                                        << (*resultHFNose)[layer].etaBin(lc.eta()) << ", "
                                        << (*resultHFNose)[layer].phiBin(lc.phi()) << "] for layer: " << layer;
    } else if (doBarrel_ && isBarrelLC) {
      resultBarrel->fill(layer, lc.eta(), lc.phi(), lcId);
      LogDebug("TICLLayerTileProducer") << "Adding layerClusterId: " << lcId << " into bin [eta,phi]: [ "
                                        << (*resultBarrel)[layer].etaBin(lc.eta()) << ", "
                                        << (*resultBarrel)[layer].phiBin(lc.phi()) << "] for layer: " << layer;
    } else if (!isBarrelLC) {
      result->fill(layer, lc.eta(), lc.phi(), lcId);
      LogDebug("TICLLayerTileProducer") << "Adding layerClusterId: " << lcId << " into bin [eta,phi]: [ "
                                        << (*result)[layer].etaBin(lc.eta()) << ", "
                                        << (*result)[layer].phiBin(lc.phi()) << "] for layer: " << layer;
    }
    lcId++;
  }
}

void TICLLayerTileProducer::produceRecHitTiles(edm::Event &evt, TICLLayerTiles* result) {
  edm::Handle<HGCRecHitCollection> ee_hits_h, fh_hits_h, bh_hits_h, hfnose_hits_h;

  auto fillFromRecHits = [this](TICLLayerTiles* result, const HGCRecHitCollection& hits) {    
    int hitId = offset_;
    for (auto const &hit : hits) {
      int layer = rhtools_.getLayerWithOffset(hit.detid()) +
      rhtools_.lastLayer(doNose_) * ((rhtools_.zside(hit.detid()) + 1) >> 1) - 1;
      float eta = rhtools_.getEta(hit.detid());
      float phi = rhtools_.getPhi(hit.detid());
      result->fill(layer, eta, phi, hitId);
      hitId++;
    }
    offset_ = hitId;
  };

  if (doNose_) {
    evt.getByToken(hgcalRecHitsHFNoseToken_, hfnose_hits_h);
    fillFromRecHits(result,*hfnose_hits_h);
  } else {
    evt.getByToken(hgcalRecHitsEEToken_, ee_hits_h);
    evt.getByToken(hgcalRecHitsFHToken_, fh_hits_h);
    evt.getByToken(hgcalRecHitsBHToken_, bh_hits_h);
    fillFromRecHits(result,*ee_hits_h);
    fillFromRecHits(result,*fh_hits_h);
    fillFromRecHits(result,*bh_hits_h);
  }
}

TICLLayerTileProducer::TICLLayerTileProducer(const edm::ParameterSet &ps)
    : detector_(ps.getParameter<std::string>("detector")), 
      isLC_(ps.getParameter<bool>("isLC")) {
  ticlGeomToken_ = esConsumes<TICLGeomHost, CaloGeometryRecord, edm::Transition::BeginRun>(edm::ESInputTag("", ""));
  ticlGeomLookupToken_ =
      esConsumes<TICLGeomLookupHost, CaloGeometryRecord, edm::Transition::BeginRun>(edm::ESInputTag("", ""));
  ticlGeomLayersToken_ =
      esConsumes<TICLGeomLayersHost, CaloGeometryRecord, edm::Transition::BeginRun>(edm::ESInputTag("", ""));
  hgcalRecHitsEEToken_ = consumes<HGCRecHitCollection>(ps.getParameter<edm::InputTag>("HGCEEInput"));
  hgcalRecHitsFHToken_ = consumes<HGCRecHitCollection>(ps.getParameter<edm::InputTag>("HGCFHInput"));
  hgcalRecHitsBHToken_ = consumes<HGCRecHitCollection>(ps.getParameter<edm::InputTag>("HGCBHInput"));  
  hgcalRecHitsHFNoseToken_ = consumes<HGCRecHitCollection>(ps.getParameter<edm::InputTag>("HGCHFNoseInput"));

  doNose_ = (detector_ == "HFNose");
  doBarrel_ = (detector_ == "Barrel");

  if (doNose_) {
    clusters_HFNose_token_ =
        consumes<std::vector<reco::CaloCluster>>(ps.getParameter<edm::InputTag>("layer_HFNose_clusters"));
    produces<TICLLayerTilesHFNose>();
  } else {
    if (doBarrel_) {
      produces<TICLLayerTilesBarrel>("ticlLayerTilesBarrel");
    }
    clusters_token_ = consumes<std::vector<reco::CaloCluster>>(ps.getParameter<edm::InputTag>("layer_clusters"));
    produces<TICLLayerTiles>();
  }
  produces<std::vector<float>>();
}

void TICLLayerTileProducer::produce(edm::Event &evt, const edm::EventSetup &) {
  std::unique_ptr<TICLLayerTilesHFNose> resultHFNose;
  std::unique_ptr<TICLLayerTiles> result;
  std::unique_ptr<TICLLayerTilesBarrel> resultBarrel;
  if (doNose_) {
    resultHFNose = std::make_unique<TICLLayerTilesHFNose>();
  } else {
    if (doBarrel_)
      resultBarrel = std::make_unique<TICLLayerTilesBarrel>();
    result = std::make_unique<TICLLayerTiles>();
  }

  offset_ = 0;

  if (isLC_) {
    produceLCTiles(evt, result.get(), resultBarrel.get(), resultHFNose.get());
  } else {
    produceRecHitTiles(evt, result.get());
  }

  if (doNose_) evt.put(std::move(resultHFNose));
  else {
    if (doBarrel_) evt.put(std::move(resultBarrel), "ticlLayerTilesBarrel");
    evt.put(std::move(result));
  }
}

void TICLLayerTileProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::string>("detector", "HGCAL");
  desc.add<edm::InputTag>("layer_clusters", edm::InputTag("hgcalMergeLayerClusters"));
  desc.add<edm::InputTag>("layer_HFNose_clusters", edm::InputTag("hgcalLayerClustersHFNose"));
  desc.add<bool>("isLC", true);
  desc.add<edm::InputTag>("HGCEEInput", edm::InputTag("HGCalRecHit", "HGCEERecHits"));
  desc.add<edm::InputTag>("HGCFHInput", edm::InputTag("HGCalRecHit", "HGCHEFRecHits"));
  desc.add<edm::InputTag>("HGCBHInput", edm::InputTag("HGCalRecHit", "HGCHEBRecHits"));
  desc.add<edm::InputTag>("HGCHFNoseInput", edm::InputTag("HGCalRecHit", "HGCHFNoseRecHits"));
  descriptions.add("ticlLayerTileProducer", desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TICLLayerTileProducer);