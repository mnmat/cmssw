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

#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticleFwd.h"

#include "SimDataFormats/CaloAnalysis/interface/SimClusterFwd.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"

struct TileParameters{
  int layer;
  double eta, phi;
};

class TICLLayerTileProducer : public edm::stream::EDProducer<> {
public:
  explicit TICLLayerTileProducer(const edm::ParameterSet &ps);
  ~TICLLayerTileProducer() override{};
  void beginRun(edm::Run const &, edm::EventSetup const &) override;
  void produce(edm::Event &, const edm::EventSetup &) override;
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  TileParameters getTileParameters(const reco::CaloCluster& lc); 
  TileParameters getTileParameters(const HGCRecHit& hit); 

  template<typename T, typename U> void fillTiles(T& results, 
      const U& objects);

  edm::EDGetTokenT<std::vector<reco::CaloCluster>> clusters_token_;
  edm::EDGetTokenT<std::vector<reco::CaloCluster>> clusters_HFNose_token_;
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsEEToken_; 
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsFHToken_;
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsBHToken_;  
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsHFNoseToken_;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geometry_token_;
  hgcal::RecHitTools rhtools_;
  std::string detector_;
  bool isLC_;
  bool doNose_;

  int offset_ = 0;
};

TileParameters TICLLayerTileProducer::getTileParameters(const reco::CaloCluster& lc){
  const auto firstHitDetId = lc.hitsAndFractions()[0].first;
  int layer = rhtools_.getLayerWithOffset(firstHitDetId) +
              rhtools_.lastLayer(doNose_) * ((rhtools_.zside(firstHitDetId) + 1) >> 1) - 1;
  float eta = lc.eta();
  float phi = lc.phi();
  return TileParameters{layer,eta,phi};
}

TileParameters TICLLayerTileProducer::getTileParameters(const HGCRecHit& hit){
  int layer = rhtools_.getLayerWithOffset(hit.detid());
  float eta = rhtools_.getEta(hit.detid());
  float phi = rhtools_.getPhi(hit.detid());
  return TileParameters{layer,eta,phi};
}

template<typename T,typename U> void TICLLayerTileProducer::fillTiles(T& results, const U& objects){
  int objId = offset_;
  for (auto const &obj : objects) {
    TileParameters par = getTileParameters(obj);
    results.fill(par.layer, par.eta, par.phi, objId);
    LogDebug("TICLLayerTileProducer") << "Adding objectId: " << objId << " into bin [eta,phi]: [ "
                                      << (results)[par.layer].etaBin(par.eta) << ", " << (results)[par.layer].phiBin(par.phi)
                                      << "] for layer: " << par.layer << std::endl;
    objId++;
  }
  offset_ = objId;
}

TICLLayerTileProducer::TICLLayerTileProducer(const edm::ParameterSet &ps)
    : detector_(ps.getParameter<std::string>("detector")), 
      isLC_(ps.getParameter<bool>("isLC")) {

  clusters_HFNose_token_ = consumes<std::vector<reco::CaloCluster>>(ps.getParameter<edm::InputTag>("layer_HFNose_clusters"));
  clusters_token_ = consumes<std::vector<reco::CaloCluster>>(ps.getParameter<edm::InputTag>("layer_clusters"));
  hgcalRecHitsEEToken_ = consumes<HGCRecHitCollection>(ps.getParameter<edm::InputTag>("HGCEEInput"));
  hgcalRecHitsFHToken_ = consumes<HGCRecHitCollection>(ps.getParameter<edm::InputTag>("HGCFHInput"));
  hgcalRecHitsBHToken_ = consumes<HGCRecHitCollection>(ps.getParameter<edm::InputTag>("HGCBHInput"));  
  hgcalRecHitsHFNoseToken_ = consumes<HGCRecHitCollection>(ps.getParameter<edm::InputTag>("HGCHFNoseInput"));
  geometry_token_ = esConsumes<CaloGeometry, CaloGeometryRecord, edm::Transition::BeginRun>();
  doNose_ = (detector_ == "HFNose");

  if (doNose_)
    produces<TICLLayerTilesHFNose>("RecHitTilesHFNose");
  else
    produces<TICLLayerTiles>("RecHitTiles");
}

void TICLLayerTileProducer::beginRun(edm::Run const &, edm::EventSetup const &es) {
  edm::ESHandle<CaloGeometry> geom = es.getHandle(geometry_token_);
  rhtools_.setGeometry(*geom);
}

void TICLLayerTileProducer::produce(edm::Event &evt, const edm::EventSetup &) {

  auto result = std::make_unique<TICLLayerTiles>();
  auto resultHFNose = std::make_unique<TICLLayerTilesHFNose>();
  offset_ = 0;
   
  edm::Handle<std::vector<reco::CaloCluster>> cluster_h;
  edm::Handle<HGCRecHitCollection> ee_hits_h;
  edm::Handle<HGCRecHitCollection> fh_hits_h;
  edm::Handle<HGCRecHitCollection> bh_hits_h;
  edm::Handle<HGCRecHitCollection> hfnose_hits_h;

  if (isLC_){
    doNose_ ? evt.getByToken(clusters_HFNose_token_, cluster_h) : evt.getByToken(clusters_token_, cluster_h);
    doNose_ ? fillTiles<TICLLayerTilesHFNose,std::vector<reco::CaloCluster>>(*resultHFNose, *cluster_h) : fillTiles<TICLLayerTiles, std::vector<reco::CaloCluster>>(*result, *cluster_h);
  } else {
    if (doNose_){
      evt.getByToken(hgcalRecHitsHFNoseToken_, hfnose_hits_h);
      fillTiles<TICLLayerTiles, HGCRecHitCollection>(*result, *hfnose_hits_h);
    }
    else{
      evt.getByToken(hgcalRecHitsEEToken_, ee_hits_h);
      evt.getByToken(hgcalRecHitsFHToken_, fh_hits_h);
      evt.getByToken(hgcalRecHitsBHToken_, bh_hits_h);

      fillTiles<TICLLayerTiles, HGCRecHitCollection>(*result, *ee_hits_h);
      fillTiles<TICLLayerTiles, HGCRecHitCollection>(*result, *fh_hits_h);
      fillTiles<TICLLayerTiles, HGCRecHitCollection>(*result, *bh_hits_h);
    }
  }

  if (doNose_){
    evt.put(std::move(resultHFNose),"RecHitTilesHFNose");
  }
  else{
    evt.put(std::move(result),"RecHitTiles");
  }
}
void TICLLayerTileProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::string>("detector", "HGCAL");
  desc.add<edm::InputTag>("layer_clusters", edm::InputTag("hgcalLayerClusters"));
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