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



class TICLLayerTileProducer : public edm::stream::EDProducer<> {
public:
  explicit TICLLayerTileProducer(const edm::ParameterSet &ps);
  ~TICLLayerTileProducer() override{};
  void beginRun(edm::Run const &, edm::EventSetup const &) override;
  void produce(edm::Event &, const edm::EventSetup &) override;
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  virtual void fillHitMap(std::map<DetId, const HGCRecHit*>& hitMap,
      const HGCRecHitCollection& recHitsEE,
      const HGCRecHitCollection& recHitsFH,
      const HGCRecHitCollection& recHitsBH) const;
  virtual void fillRecHits(std::vector<const HGCRecHit*>& recHits,
      const HGCRecHitCollection& recHitsEE,
      const HGCRecHitCollection& recHitsFH,
      const HGCRecHitCollection& recHitsBH) const;
  //template <class T, class U>
  //void fillTiles(T& objects, U& results);
  edm::EDGetTokenT<std::vector<reco::CaloCluster>> clusters_token_;
  edm::EDGetTokenT<std::vector<reco::CaloCluster>> clusters_HFNose_token_;
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsEEToken_; 
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsFHToken_;
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsBHToken_;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geometry_token_;
  //edm::EDGetTokenT<std::vector<CaloParticle>> caloParticlesToken_;
  hgcal::RecHitTools rhtools_;
  std::string detector_;
  bool isLC_;
  bool doNose_;
};

void TICLLayerTileProducer::fillRecHits(std::vector<const HGCRecHit*>& recHits,
                                const HGCRecHitCollection& rechitsEE,
                                const HGCRecHitCollection& rechitsFH,
                                const HGCRecHitCollection& rechitsBH) const {
  recHits.clear();
  for (const auto& hit : rechitsEE) {
    recHits.push_back(&hit);
  }

  for (const auto& hit : rechitsFH) {
    recHits.push_back(&hit);
  }

  for (const auto& hit : rechitsBH) {
    recHits.push_back(&hit);
  }
} // end of EfficiencyStudies::fillHitMap

void TICLLayerTileProducer::fillHitMap(std::map<DetId,const HGCRecHit*>& hitMap,
                                const HGCRecHitCollection& rechitsEE,
                                const HGCRecHitCollection& rechitsFH,
                                const HGCRecHitCollection& rechitsBH) const {
  hitMap.clear();
  for (const auto& hit : rechitsEE) {
    hitMap.emplace(hit.detid(), &hit);
  }

  for (const auto& hit : rechitsFH) {
    hitMap.emplace(hit.detid(), &hit);
  }

  for (const auto& hit : rechitsBH) {
    hitMap.emplace(hit.detid(), &hit);
  }
} // end of EfficiencyStudies::fillHitMap

/*
template <class T, class U>
void TICLLayerTileProducer::fillTiles(T& objects, U& result){
  int objId = 0;
  int layer;
  float eta, phi;
  for (auto const &obj : objects) {
    if (typeid(obj)){
      const auto firstHitDetId = obj.hitsAndFractions()[0].first;
      layer = rhtools_.getLayerWithOffset(firstHitDetId) +
                  rhtools_.lastLayer(doNose_) * ((rhtools_.zside(firstHitDetId) + 1) >> 1) - 1;
      eta = obj.eta();
      phi = obj.phi();
    }
    else {
      layer = rhtools_.getLayerWithOffset(obj->detid());
      eta = rhtools_.getEta(obj->detid());
      phi = rhtools_.getPhi(obj->detid());
    }


    assert(layer >= 0);
    result->fill(layer, eta, phi, objId);
    LogDebug("TICLLayerTileProducer") << "Adding layerClusterId: " << objId << " into bin [eta,phi]: [ "
                                      << (*result)[layer].etaBin(eta) << ", " << (*result)[layer].phiBin(phi)
                                      << "] for layer: " << layer << std::endl;
    objId++;
  }
}
*/
TICLLayerTileProducer::TICLLayerTileProducer(const edm::ParameterSet &ps)
    : detector_(ps.getParameter<std::string>("detector")), 
      isLC_(ps.getParameter<bool>("isLC")) {
  clusters_HFNose_token_ =
      consumes<std::vector<reco::CaloCluster>>(ps.getParameter<edm::InputTag>("layer_HFNose_clusters"));
  clusters_token_ = consumes<std::vector<reco::CaloCluster>>(ps.getParameter<edm::InputTag>("layer_clusters"));
  hgcalRecHitsEEToken_ = consumes<HGCRecHitCollection>(ps.getParameter<edm::InputTag>("HGCEEInput"));
  hgcalRecHitsFHToken_ = consumes<HGCRecHitCollection>(ps.getParameter<edm::InputTag>("HGCFHInput"));
  hgcalRecHitsBHToken_ = consumes<HGCRecHitCollection>(ps.getParameter<edm::InputTag>("HGCBHInput"));
  //caloParticlesToken_ = consumes<std::vector<CaloParticle>>(ps.getParameter<edm::InputTag>("caloParticles")); 
  geometry_token_ = esConsumes<CaloGeometry, CaloGeometryRecord, edm::Transition::BeginRun>();

  doNose_ = (detector_ == "HFNose");

  if (doNose_)
    produces<TICLLayerTilesHFNose>("Test");
  else
    produces<TICLLayerTiles>("Test");
  produces<std::vector<float>>("Test");
}

void TICLLayerTileProducer::beginRun(edm::Run const &, edm::EventSetup const &es) {
  edm::ESHandle<CaloGeometry> geom = es.getHandle(geometry_token_);
  rhtools_.setGeometry(*geom);
}

void TICLLayerTileProducer::produce(edm::Event &evt, const edm::EventSetup &) {

  std::cout << "in TICLLayerTileProducer.cc, start putting RecHits in event" << std::endl;

  auto result = std::make_unique<TICLLayerTiles>();
  auto resultHFNose = std::make_unique<TICLLayerTilesHFNose>();
  auto test = std::make_unique<std::vector<float>>();
  
  edm::Handle<std::vector<reco::CaloCluster>> cluster_h;
  edm::Handle<HGCRecHitCollection> ee_hits;
  edm::Handle<HGCRecHitCollection> fh_hits;
  edm::Handle<HGCRecHitCollection> bh_hits;
  std::map<DetId, const HGCRecHit*> hitMap;
  std::vector<const HGCRecHit*> recHits;

  /*
  edm::Handle<std::vector<CaloParticle>> CaloParticles;
  evt.getByToken(caloParticlesToken_, CaloParticles);
  const CaloParticleCollection& cps = *CaloParticles;

  */
  /*
  if (doNose_)
    evt.getByToken(clusters_HFNose_token_, cluster_h);
  else
    evt.getByToken(clusters_token_, cluster_h);
  */

  if (isLC_){
    doNose_ ? evt.getByToken(clusters_HFNose_token_, cluster_h) : evt.getByToken(clusters_token_, cluster_h);
    //const auto &objects = *cluster_h;
    //doNose_ ? fillTiles(*cluster_h, result) : fillTiles(*cluster_h, resultHFNose);

    int lcId = 0;
    for (auto const &lc : *cluster_h) {
      const auto firstHitDetId = lc.hitsAndFractions()[0].first;
      int layer = rhtools_.getLayerWithOffset(firstHitDetId) +
                  rhtools_.lastLayer(doNose_) * ((rhtools_.zside(firstHitDetId) + 1) >> 1) - 1;
      assert(layer >= 0);

      doNose_ ? resultHFNose->fill(layer, lc.eta(), lc.phi(), lcId) : result->fill(layer, lc.eta(), lc.phi(), lcId);

      LogDebug("TICLLayerTileProducer") << "Adding layerClusterId: " << lcId << " into bin [eta,phi]: [ "
                                        << (*result)[layer].etaBin(lc.eta()) << ", " << (*result)[layer].phiBin(lc.phi())
                                        << "] for layer: " << layer << std::endl;
      lcId++;
    }
  }
  else {
    evt.getByToken(hgcalRecHitsEEToken_, ee_hits);
    evt.getByToken(hgcalRecHitsFHToken_, fh_hits);
    evt.getByToken(hgcalRecHitsBHToken_, bh_hits);
    fillHitMap(hitMap, *ee_hits, *fh_hits, *bh_hits);
    fillRecHits(recHits, *ee_hits, *fh_hits, *bh_hits);
    //doNose_ ? fillTiles(recHits, result) : fillTiles(recHits, resultHFNose);

    //doNose_ ? evt.getByToken(clusters_HFNose_token_, cluster_h) : evt.getByToken(clusters_token_, cluster_h);
    //const auto &objects = *cluster_h;
    //doNose_ ? fillTiles(*cluster_h, result) : fillTiles(*cluster_h, resultHFNose);

    int objId = 0;
    for (auto const &rhit : recHits) {
      int layer = rhtools_.getLayerWithOffset(rhit->detid());
      float eta = rhtools_.getEta(rhit->detid());
      float phi = rhtools_.getPhi(rhit->detid());

      assert(layer >= 0);

      objId = rhit->detid().rawId();

      if (doNose_){
        resultHFNose->fill(layer, eta, phi, objId);
      }
      else {
        result->fill(layer, eta, phi, objId);
      }
      LogDebug("TICLLayerTileProducer") << "Adding layerClusterId: " << objId << " into bin [eta,phi]: [ "
                                        << (*result)[layer].etaBin(eta) << ", " << (*result)[layer].phiBin(phi)
                                        << "] for layer: " << layer << std::endl;
      objId++;
    }
    /*
    int objId=0;
    for (const auto& it_cp : cps ) {
      const CaloParticle& cp = ((it_cp));
      const SimClusterRefVector&  simclusters = cp.simClusters(); 
      for (const auto& it_simc : simclusters) {
        const SimCluster& simc = (*(it_simc));
        const auto& sc_haf = simc.hits_and_fractions();
        for (const auto&  it_sc_haf : sc_haf) {
          DetId detid_ = (it_sc_haf.first);
          std::map<DetId,const HGCRecHit *>::const_iterator itcheck = hitMap.find(detid_);
          if (itcheck != hitMap.end()) {
            int layer = rhtools_.getLayerWithOffset(detid_);
            float eta = rhtools_.getEta(detid_);
            float phi = rhtools_.getPhi(detid_);

            assert(layer >= 0);

            if (doNose_){
              resultHFNose->fill(layer, eta, phi, objId);
            }
            else {
              result->fill(layer, eta, phi, objId);
            }
            LogDebug("TICLLayerTileProducer") << "Adding recHitId: " << objId << " into bin [eta,phi]: [ "
                                              << (*result)[layer].etaBin(eta) << ", " << (*result)[layer].phiBin(phi)
                                              << "] for layer: " << layer << std::endl;
            objId++;
          }
        }
      }
    }
    */
  }
  
  //const auto &object = (isLC_? recHits : *cluster_h)
  /*

  int objId = 0;
  for (auto const &obj : objects) {
    if (isLC_){
      const auto firstHitDetId = lc.hitsAndFractions()[0].first;
      int layer = rhtools_.getLayerWithOffset(firstHitDetId) +
                  rhtools_.lastLayer(doNose_) * ((rhtools_.zside(firstHitDetId) + 1) >> 1) - 1;
    }
    else {
      int layer = obj.layer();
    }

    assert(layer >= 0);

    if (doNose_){
      resultHFNose->fill(layer, obj.eta(), obj.phi(), objId);
    }
    else {
      result->fill(layer, obj.eta(), obj.phi(), objId);
    }
    LogDebug("TICLLayerTileProducer") << "Adding layerClusterId: " << objId << " into bin [eta,phi]: [ "
                                      << (*result)[layer].etaBin(obj.eta()) << ", " << (*result)[layer].phiBin(obj.phi())
                                      << "] for layer: " << layer << std::endl;
    objId++;
  }
  */

  std::cout << "put" << std::endl;
  if (doNose_){
    evt.put(std::move(resultHFNose),"Test");
  }
  else{
    evt.put(std::move(result),"Test");
  }
  evt.put(std::move(test),"Test");
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
  //desc.add<edm::InputTag>("caloParticles", edm::InputTag("mix", "MergedCaloTruth"));
  descriptions.add("ticlLayerTileProducer", desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TICLLayerTileProducer);
