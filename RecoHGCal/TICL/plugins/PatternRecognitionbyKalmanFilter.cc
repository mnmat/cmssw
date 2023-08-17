// Author: Mark Matthewman - mark.matthewman@cern.ch
// Date: 05/2022
#include <algorithm>
#include <set>
#include <vector>
#include <typeinfo>

#include "PatternRecognitionbyKalmanFilter.h"

#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackingRecHit/interface/InvalidTrackingRecHit.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/HGCalReco/interface/HGCTrackingRecHit.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimatorBase.h"
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h"
#include "TrackingTools/KalmanUpdators/interface/KFUpdator.h"
#include "TrackingTools/PatternTools/interface/TrajectoryStateUpdator.h"
#include "TrackingTools/PatternTools/interface/TrajMeasLessEstim.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

#include "HGCTracker.h"

using namespace ticl;

template <typename TILES>
PatternRecognitionbyKalmanFilter<TILES>::PatternRecognitionbyKalmanFilter(const edm::ParameterSet &conf, edm::ConsumesCollector iC)
    : PatternRecognitionAlgoBaseT<TILES>(conf, iC),
      caloGeomToken_(iC.esConsumes<CaloGeometry, CaloGeometryRecord>()),
      radlen_(conf.getParameter<std::vector<double>>("radlen")),
      xi_(conf.getParameter<std::vector<double>>("xi")),
      propName_(conf.getParameter<std::string>("propagator")),
      propNameOppo_(conf.getParameter<std::string>("propagatorOpposite")),
      bfieldtoken_(iC.esConsumes<MagneticField, IdealMagneticFieldRecord>()),
      propagatortoken_(iC.esConsumes<Propagator, TrackingComponentsRecord>(edm::ESInputTag("",propName_))),
      propagatorOppoToken_(iC.esConsumes<Propagator, TrackingComponentsRecord>(edm::ESInputTag("",propNameOppo_))),
      estimatorToken_(iC.esConsumes<Chi2MeasurementEstimatorBase, TrackingComponentsRecord>(edm::ESInputTag("","Chi2"))),
      updatorToken_(iC.esConsumes<TrajectoryStateUpdator, TrackingComponentsRecord>(edm::ESInputTag("","KFUpdator"))),
      trackToken_(iC.consumes<reco::TrackCollection>(conf.getParameter<edm::InputTag>("tracks"))),
      hgcalRecHitsEEToken_(iC.consumes<HGCRecHitCollection>(conf.getParameter<edm::InputTag>("HGCEEInput"))),
      hgcalRecHitsFHToken_(iC.consumes<HGCRecHitCollection>(conf.getParameter<edm::InputTag>("HGCFHInput"))),
      hgcalRecHitsBHToken_(iC.consumes<HGCRecHitCollection>(conf.getParameter<edm::InputTag>("HGCBHInput"))),
      diskToken_(iC.esConsumes<HGCDiskGeomDetVector, CaloGeometryRecord>()),
      hgcTrackerToken_(iC.esConsumes<HGCTracker, CaloGeometryRecord>()),
      rescaleFTSError_(conf.getParameter<double>("rescaleFTSError")),
      scaleWindow_(conf.getParameter<double>("scaleWindow")),
      geomCacheId_(0)
      {};


template <typename TILES>
std::pair<float,float> PatternRecognitionbyKalmanFilter<TILES>::covarianceTransform(const TrajectoryStateOnSurface &tsos){
  // Get Position
  double x = tsos.globalPosition().x();
  double y = tsos.globalPosition().y();
  double z = tsos.globalPosition().z();
  
  // Calculate Jacobian
  AlgebraicMatrix22 theJacobian;
  double sqrt_term = std::sqrt((x*x + y*y) / (z*z) + 1);
  double denom_eta = (x*x + y*y) * (x*x + y*y + z*z);
  theJacobian(0,0) = - (x*z*z * sqrt_term) / denom_eta; // deta_dx
  theJacobian(0,1) = - (y*z*z * sqrt_term) / denom_eta; // deta_dy

  double denom_phi = x*x + y*y;
  theJacobian(1,0) = - y / denom_phi; // dphi_dx
  theJacobian(1,1) = x / denom_phi; // dphi_dy

  // Get covariance in x-y coordinates
  AlgebraicMatrix55 localErrorMatrix = tsos.localError().matrix();
  AlgebraicMatrix22 covMatrixXY;
  covMatrixXY(0,0) = localErrorMatrix(3,3);
  covMatrixXY(0,1) = localErrorMatrix(3,4);
  covMatrixXY(1,0) = localErrorMatrix(3,4);
  covMatrixXY(1,1) = localErrorMatrix(4,4);

  // Transform LocalError from x-y coordinates to eta-phi
  // TODO: Look for more suitable functions in ROOT::Math namespace
  AlgebraicMatrix22 covMatrixEtaPhi = ROOT::Math::Transpose(theJacobian) * covMatrixXY * theJacobian;
  
  return std::pair<float,float>{covMatrixEtaPhi(0,0),covMatrixEtaPhi(1,1)};
}


template<typename TILES>
template<class Start> 
std::vector<TempTrajectory>
PatternRecognitionbyKalmanFilter<TILES>::advanceOneLayer(const Start &start, 
                                              const HGCDiskLayer * diskLayer, 
                                              const TILES &tiles, PropagationDirection direction, 
                                              bool &isSilicon, 
                                              TempTrajectory traj){  
                                                                            
  auto disk = isSilicon? diskLayer->second : diskLayer->first;         
  std::vector<TempTrajectory> ret;
  const Propagator &prop = (direction == alongMomentum ? *propagator_ : *propagatorOppo_);
  TrajectoryStateOnSurface tsos = prop.propagate(start, disk->surface());
  if (!tsos.isValid()) return ret; 

  // Collect hits with estimate
  int layer = disk->layer();
  std::vector<std::shared_ptr<HGCTrackingRecHit>> hitptrs;

  if (disk->zside() > 0){
    auto lastLayerPerSide = static_cast<int>(rhtools_.lastLayer(false));
    hitptrs = measurements(tsos, tiles, layer+lastLayerPerSide-1);
  }
  else { 
    hitptrs = measurements(tsos, tiles, layer-1);
  }
  // Sort according to isSilicon
  /*
  auto ptrSortSilicon = [](const HGCTrackingRecHit a, const HGCTrackingRecHit b) -> bool { return a.isSilicon < b.isSilicon; };
  auto ptrSortScintillator = [](const HGCTrackingRecHit a, const HGCTrackingRecHit b) -> bool { return a.isSilicon > b.isSilicon; };

  if (isSilicon){
    std::sort(hitptrs.begin(), hitptrs.end(),ptrSortSilicon);
  } else {
    std::sort(hitptrs.begin(), hitptrs.end(),ptrSortScintillator);
  }
  */

  // Create Trajectory Measurements for each Measurement with the correct TSOS
  std::vector<TrajectoryMeasurement> meas;
  for(const auto &hitptr : hitptrs){
    if (hitptr->isValid()){
      auto id = static_cast<int32_t>(hitptr->rawId());
      if (rhtools_.isSilicon(id) != isSilicon){
        isSilicon = !isSilicon;
        disk = isSilicon? diskLayer->second : diskLayer->first;
        if (!disk){
          std::cout << "Disk doesn't exist" << std::endl;
        }
        tsos = prop.propagate(start, disk->surface());
      }
    }
    auto testid = static_cast<int32_t>(hitptr->rawId());
    auto mest_pair = (*estimator_).estimate(tsos,*hitptr);
    if(mest_pair.first){
      meas.emplace_back(tsos,hitptr,mest_pair.second);
    }
  }

  // Sort according to chi2
  std::sort(meas.begin(), meas.end(),TrajMeasLessEstim());

  for (const TrajectoryMeasurement &tm : meas){
    TrajectoryStateOnSurface updated = updator_->update(tm.forwardPredictedState(),*tm.recHit());
    ret.push_back(traj.foundHits() ? traj: TempTrajectory(traj.direction(),0));
    ret.back().push(TrajectoryMeasurement(tm.forwardPredictedState(),
                                          updated,
                                          tm.recHit(),
                                          tm.estimate()),
                    tm.estimate());
  }

  if (!meas.empty()) return ret;
  auto missing = TrackingRecHit::missing;
  ret.push_back(traj.foundHits()? traj : TempTrajectory(traj.direction(),0));
  ret.back().push(TrajectoryMeasurement(tsos, std::make_shared<InvalidTrackingRecHit>(*disk,missing)));

  return ret;
}

template <typename TILES>
void PatternRecognitionbyKalmanFilter<TILES>::mergeRecHitCollections(std::vector<const HGCRecHit*>& recHitCollection,
                                const HGCRecHitCollection& rechitsEE,
                                const HGCRecHitCollection& rechitsFH,
                                const HGCRecHitCollection& rechitsBH) const {
  
  recHitCollection.clear();
  for (const auto& hit : rechitsEE) {
    recHitCollection.push_back(&hit);
  }

  for (const auto& hit : rechitsFH) {
    recHitCollection.push_back(&hit);
  }

  for (const auto& hit : rechitsBH) {
    recHitCollection.push_back(&hit);
  }
}

template <typename TILES>
std::vector<std::shared_ptr<HGCTrackingRecHit>> PatternRecognitionbyKalmanFilter<TILES>::measurements(
      const TrajectoryStateOnSurface &tsos, 
      const TILES &tiles, 
      int layer){
  
  std::vector<std::shared_ptr<HGCTrackingRecHit>> ret;

  // define search window and get bins

  float eta = tsos.globalPosition().eta();
  float phi = tsos.globalPosition().phi();

  std::pair<float,float> localErrorEtaPhi  = covarianceTransform(tsos);

  float etaMin = eta - scaleWindow_*std::max(etaBinSize,localErrorEtaPhi.first);
  float etaMax = eta + scaleWindow_*std::max(etaBinSize,localErrorEtaPhi.first);
  float phiMin = phi - scaleWindow_*std::max(phiBinSize,localErrorEtaPhi.second);
  float phiMax = phi + scaleWindow_*std::max(phiBinSize,localErrorEtaPhi.second);

  auto bins = tiles[layer].searchBoxEtaPhi(etaMin, etaMax, phiMin, phiMax);


  // loop over candidates

  for (int ieta = bins[0]; ieta <= bins[1]; ieta++) {
    auto offset = ieta * nPhiBin;
    for (int phi = bins[2]; phi <= bins[3]; phi++) {
      int iphi = ((phi % nPhiBin + nPhiBin) % nPhiBin);
      if (!tiles[layer][offset + iphi].empty()) {
        for(auto hit: tiles[layer][offset + iphi]) {
          const auto rec = *recHitCollection[hit];
          const auto detid = rec.detid();

          GlobalPoint globalpoint = rhtools_.getPosition(detid);
          LocalPoint localpoint = tsos.surface().toLocal(globalpoint);
          bool isSilicon = rhtools_.isSilicon(detid); 

          auto hitptr = std::make_shared<HGCTrackingRecHit>(detid,localpoint,rhtools_.getLocalError(detid));
          ret.push_back(hitptr);
        }
      }
    }
  }
return ret;
}

template <typename TILES>
void PatternRecognitionbyKalmanFilter<TILES>::init(
    const edm::Event& ev, const edm::EventSetup& es){

    //Get Calo Geometry
    if (es.get<CaloGeometryRecord>().cacheIdentifier() != geomCacheId_) {
      geomCacheId_ = es.get<CaloGeometryRecord>().cacheIdentifier();
      const CaloGeometry* geom = &es.getData(caloGeomToken_);
      rhtools_.setGeometry(*geom);
      hgcTracker_ = &es.getData(hgcTrackerToken_);
    } 

    bfield_ = es.getHandle(bfieldtoken_);
    propagator_ = es.getHandle(propagatortoken_); // Can I move this up?
    propagatorOppo_ = es.getHandle(propagatorOppoToken_); // Can I move this up?
    estimator_ = es.getHandle(estimatorToken_); // Can I move this up?
    updator_ = es.getHandle(updatorToken_); // Can I move this up?

    edm::Handle<HGCRecHitCollection> ee_hits;
    edm::Handle<HGCRecHitCollection> fh_hits;
    edm::Handle<HGCRecHitCollection> bh_hits;

    ev.getByToken(hgcalRecHitsEEToken_, ee_hits);
    ev.getByToken(hgcalRecHitsFHToken_, fh_hits);
    ev.getByToken(hgcalRecHitsBHToken_, bh_hits);
    mergeRecHitCollections(recHitCollection, *ee_hits, *fh_hits, *bh_hits);
}


template <typename TILES>
void PatternRecognitionbyKalmanFilter<TILES>::makeTrajectories(
    const typename PatternRecognitionAlgoBaseT<TILES>::Inputs &input,
    std::vector<KFHit>& kfhits) {


  edm::EventSetup const &es = input.es;
  edm::Event const &ev = input.ev;
  init(ev,es);

  const TILES &tiles = input.tiles;
  //dumpTiles(tiles);

  // Build tsos from TrackCollection
  edm::Handle<reco::TrackCollection> tracks_h;
  ev.getByToken(trackToken_,tracks_h);

  const reco::TrackCollection& tkx = *tracks_h; 
  if (tkx.empty()){
    edm::LogWarning("PatternRecognitionbyKalmanFilter") << "No track to extrapolate from first disk found! Exited PatternRecognitionbyKalmanFilter" << std::endl;
    return;
  }

  trackId = 0;
  for(auto tk: tkx){
    /*
    if (std::abs(tk.p())<97){
      trackId++;
      continue;
    }
    */

    FreeTrajectoryState fts = trajectoryStateTransform::outerFreeState(tk,bfield_.product());
    if (rescaleFTSError_!=1){
      //TODO: Message Logging for Error
      fts.rescaleError(rescaleFTSError_);
    }

    // Extrapolate muon position from the tracker

    int zside = fts.momentum().eta() > 0 ? +1 : -1;
    PropagationDirection direction = alongMomentum;

    const HGCDiskLayer* layerdisk = hgcTracker_->firstDisk(zside,direction);
    bool isSilicon = true;
    
    std::vector<TempTrajectory> traj = advanceOneLayer(fts, layerdisk,tiles, direction, isSilicon, TempTrajectory(direction,0));
    if (traj.empty()){
      edm::LogWarning("PatternRecognitionbyKalmanFilter") << "No valid Trajectory found! Skip track!" << std::endl; 
      trackId++;
      continue;
    }

    auto lm = traj.back().lastMeasurement();
    TrajectoryStateOnSurface tsos_kf = lm.updatedState();
    int layer = layerdisk->second->layer();
    KFHit *kfhit = new KFHit(tsos_kf, lm.recHit()->geographicalId(), tk, trackId, layer);
    kfhits.push_back(*kfhit);

    std::vector<TempTrajectory> traj_kf;
    traj_kf.push_back(traj.back());

    // Loop over all disks
    for(layerdisk = hgcTracker_->nextDisk(layerdisk,direction,isSilicon); layerdisk != nullptr; layerdisk = hgcTracker_->nextDisk(layerdisk,direction,isSilicon)){
      std::vector<TempTrajectory> newcands_kf;
      for(TempTrajectory & cand : traj_kf){
        TrajectoryStateOnSurface start = cand.lastMeasurement().updatedState();


        std::vector<TempTrajectory> hisTrajs = advanceOneLayer(start, layerdisk, tiles, direction, isSilicon, cand);
        for(TempTrajectory & t : hisTrajs){
          auto lm = t.lastMeasurement();
          newcands_kf.push_back(t);

          TrajectoryStateOnSurface tsos_kf = lm.updatedState();
          layer = layerdisk->second->layer();
          KFHit *kfhit = new KFHit(tsos_kf, lm.recHit()->geographicalId(), tk, trackId, layer);
          kfhits.push_back(*kfhit);
          break;
        }
      }
      traj_kf.swap(newcands_kf);
    }
    
    trackId++;
  }
}

template <typename TILES>
void PatternRecognitionbyKalmanFilter<TILES>::dumpTiles(const TILES &tiles) const {
  edm::LogInfo("PatternRecognitionbyKalmanFilter") << "Entered dumpTiles" << std::endl;
  constexpr int nEtaBin = TILES::constants_type_t::nEtaBins;
  constexpr int nPhiBin = TILES::constants_type_t::nPhiBins;
  edm::LogInfo("PatternRecognitionbyKalmanFilter") << nEtaBin << "\t" <<nPhiBin << std::endl;
  auto lastLayerPerSide = static_cast<int>(rhtools_.lastLayer(false));
  edm::LogInfo("PatternRecognitionbyKalmanFilter") << lastLayerPerSide << std::endl;
  int maxLayer = 2 * lastLayerPerSide - 1;
  int count = 0;
  for (int layer = 0; layer <= maxLayer; layer++) {
    for (int ieta = 0; ieta < nEtaBin; ieta++) {
      auto offset = ieta * nPhiBin;
      for (int phi = 0; phi < nPhiBin; phi++) {
        int iphi = ((phi % nPhiBin + nPhiBin) % nPhiBin);
        if (!tiles[layer][offset + iphi].empty()) {
          edm::LogInfo("PatternRecognitionbyKalmanFilter") << "Layer: " << layer << " ieta: " << ieta << " phi: " << phi
                                                         << " " << tiles[layer][offset + iphi].size() << std::endl;
          for(auto hit: tiles[layer][offset + iphi]) {
            const auto rec = *recHitCollection[hit];
            const auto detid = rec.detid();
            std::cout << "Layer: " << layer << " ieta: " << ieta << " phi: " << phi
                                                         << " DetId: " << static_cast<int>(detid.rawId()) << std::endl;          
          }
          count++;
        }
      }
    }
  }
  edm::LogInfo("PatternRecognitionbyKalmanFilter") << "Number of RecHits: " << count << std::endl;
}

template <typename TILES>
void PatternRecognitionbyKalmanFilter<TILES>::fillPSetDescription(edm::ParameterSetDescription &iDesc) {
  iDesc.add<int>("algo_verbosity", 0);
  iDesc.add<std::vector<double>>("radlen",{});
  iDesc.add<std::vector<double>>("xi",{});
  iDesc.add<std::string>("propagator", "PropagatorWithMaterial"); // Analytical Propagator 
  iDesc.add<std::string>("propagatorOpposite", "PropagatorWithMaterialOpposite");
  iDesc.add<std::string>("estimator", "Chi2");
  iDesc.add<edm::InputTag>("tracks", edm::InputTag("generalTracks"));
  iDesc.add<edm::InputTag>("HGCEEInput", edm::InputTag("HGCalRecHit", "HGCEERecHits"));
  iDesc.add<edm::InputTag>("HGCFHInput", edm::InputTag("HGCalRecHit", "HGCHEFRecHits"));
  iDesc.add<edm::InputTag>("HGCBHInput", edm::InputTag("HGCalRecHit", "HGCHEBRecHits"));
  iDesc.add<double>("rescaleFTSError",1);
  iDesc.add<double>("scaleWindow",1);
}

template class ticl::PatternRecognitionbyKalmanFilter<TICLLayerTiles>;
template class ticl::PatternRecognitionbyKalmanFilter<TICLLayerTilesHFNose>;
