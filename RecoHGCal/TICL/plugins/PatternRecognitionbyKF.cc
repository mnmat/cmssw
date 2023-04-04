// Author: Mark Matthewman - mark.matthewman@cern.ch
// Date: 05/2022
#include <algorithm>
#include <set>
#include <vector>
#include <typeinfo>

#include "PatternRecognitionbyKF.h"

#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackingRecHit/interface/InvalidTrackingRecHit.h"
#include "DataFormats/HGCTrackingRecHit/interface/HGCTrackingRecHit.h"

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

using namespace ticl;

template <typename TILES>
PatternRecognitionbyKF<TILES>::PatternRecognitionbyKF(const edm::ParameterSet &conf, edm::ConsumesCollector iC)
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
      geomCacheId_(0){
};

template<typename TILES>
void PatternRecognitionbyKF<TILES>::calculateLocalError(DetId id, const HGCalGeometry* hgcalgeom_){
  // Calculations based on normalized second moment of area
  if(rhtools_.isSilicon(id)){
    double A = hgcalgeom_->getArea(id);
    double a = sqrt(2*A/(3*sqrt(3))); // side length hexagon
    double var = pow(a,4)*5*sqrt(3)/(16*A);
    lerr[id] = LocalError(var, 0, var);
  }
  else{
    const HGCalDDDConstants* ddd = &(hgcalgeom_->topology().dddConstants());
    // Get outer and inner radius
    const GlobalPoint &pos = rhtools_.getPosition(id);
    double r = sqrt(pos.x()*pos.x() + pos.y()*pos.y());
    auto radiusLayer = ddd->getRadiusLayer(rhtools_.getLayer(id));
    int idx = static_cast<int>(std::lower_bound(radiusLayer.begin(), radiusLayer.end(),r)-radiusLayer.begin());
    double rmax = radiusLayer[idx];
    double rmin = radiusLayer[idx-1];
    // Get angles
    double phi = rhtools_.getPhi(id) + M_PI; // set to radians [0, 2pi]
    double dphi = rhtools_.getScintDEtaDPhi(id).second;
    double phimin = phi - 0.5*dphi;
    double phimax = phi + 0.5*dphi;

    // FIXME!!! Using getArea() function results in lower efficiency due to incorrect position of KF hits. To be understood and getArea() to be implemented
    // double A = hgcalgeom_->getArea(id); TOD
    double A = (rmax*rmax - rmin*rmin)*M_PI*dphi/(2*M_PI); 
    // Calculate local error
    double ex2 = 1/(8*A) * (pow(rmax,4) - pow(rmin,4)) * (-phimin - sin(phimin)*cos(phimin) + phimax + sin(phimax)*cos(phimax));
    double ex = 1/(3*A) * (pow(rmax,3) - pow(rmin,3)) * (sin(phimax) - sin(phimin));
    double varx = ex2 - ex*ex;
    double ey2 = 1/(8*A) * (pow(rmax,4) - pow(rmin,4)) * (-phimin + sin(phimin)*cos(phimin) + phimax - sin(phimax)*cos(phimax));
    double ey = 1/(3*A) * (pow(rmax,3) - pow(rmin,3)) * (cos(phimin) - cos(phimax));
    double vary = ey2 - ey*ey;
    double varxy = 1/(16*A)*(pow(rmax,4)-pow(rmin,4))*(cos(2*phimin)-cos(2*phimax)) - ex*ey;
    lerr[id] = LocalError(varx, varxy, vary);
  }
} 

template<typename TILES>
const HGCDiskGeomDet * PatternRecognitionbyKF<TILES>::switchDisk(const HGCDiskGeomDet * from, 
                                                                const std::vector<HGCDiskGeomDet *> &vec, 
                                                                bool isSilicon) const{                                                           
  auto it = std::find(vec.begin(), vec.end(),from);
  if (it == vec.end()) throw cms::Exception("LogicError", "nextDisk called with invalid starting disk");
  isSilicon? --it: ++it;
  return *(it);
}

template<typename TILES>
const HGCDiskGeomDet * PatternRecognitionbyKF<TILES>::nextDisk(const HGCDiskGeomDet * from, 
                                                              PropagationDirection direction, 
                                                              const std::vector<HGCDiskGeomDet *> &vec, 
                                                              bool isSilicon) const{
  auto it = std::find(vec.begin(), vec.end(),from);
  if (it == vec.end()) throw cms::Exception("LogicError", "nextDisk called with invalid starting disk");
  int currentLayer = (*it)->layer();
  if (direction == alongMomentum){
    if ((*it == vec.back()) || (*it == vec.rbegin()[1])) return nullptr; // if disk last silicon OR last scintillator disk
    while ((*it)->layer()<currentLayer+2){
      if ((*it) == vec.back()) break;
      ++it;
      if(((*(it))->isSilicon() == isSilicon) && ((*(it))->layer()==currentLayer+1)) return *(it);
     }
  } else{
    if (it == vec.begin()) return nullptr;
    while ((*it)->layer()>currentLayer-2){
      if ((it) == vec.begin()) break;
      --it;
      if(((*(it))->isSilicon() == isSilicon) && ((*(it))->layer()==currentLayer-1)) return *(it);
    }
  }
  return nullptr; // Returns nullptr if nextDisk reached end of detector
}

template<typename TILES>
template<class Start> 
std::vector<TempTrajectory>
PatternRecognitionbyKF<TILES>::advanceOneLayer(const Start &start, 
                                              const HGCDiskGeomDet * disk, 
                                              const std::vector<HGCDiskGeomDet *> &disks, 
                                              const TILES &tiles, PropagationDirection direction, 
                                              bool &isSilicon, 
                                              TempTrajectory traj){

  std::vector<TempTrajectory> ret;

  const Propagator &prop = (direction == alongMomentum ? *propagator_ : *propagatorOppo_);
  TrajectoryStateOnSurface tsos = prop.propagate(start, disk->surface());
  if (!tsos.isValid()) return ret;
  float r = sqrt(pow(tsos.globalPosition().x(),2)+pow(tsos.globalPosition().y(),2));
  if (((disk->rmin() > r) && (!isSilicon)) || (((r > disk->rmax()) && (isSilicon)))) {
    disk = switchDisk(disk, disks, isSilicon);
    isSilicon = !isSilicon;
    tsos = prop.propagate(start, disk->surface());
  }

  // Collect hits with estimate
  int depth = disk->layer()+1;
  auto meas = measurements(tsos, *estimator_, tiles, depth);
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

template<typename TILES>
void PatternRecognitionbyKF<TILES>::makeDisks(int subdet, int disks, const CaloGeometry* geom_) {

    const CaloSubdetectorGeometry *subGeom = geom_->getSubdetectorGeometry(DetId::Detector(subdet), ForwardSubdetector::ForwardEmpty);
    auto geomEE = static_cast<const HGCalGeometry*>(subGeom);
    const HGCalDDDConstants* ddd = &(geomEE->topology().dddConstants());

    std::vector<float>  rmax(disks, 0), rmin(disks, 9e9);
    std::vector<double> zsumPos(disks), zsumNeg(disks);
    std::vector<int> countPos(disks), countNeg(disks);
    const std::vector<DetId> & ids = subGeom->getValidDetIds();
  
    for (auto & i : ids) {
        calculateLocalError(i,geomEE);

        const GlobalPoint & pos = rhtools_.getPosition(i); 
        int layer = rhtools_.getLayer(i)-1;
        float z = pos.z();
        float rho = pos.perp();
        int side = z > 0 ? +1 : -1;

        (side > 0 ? zsumPos : zsumNeg)[layer] += z;
        (side > 0 ? countPos : countNeg)[layer]++;
        if (rho > rmax[layer]) rmax[layer] = rho;
        if (rho < rmin[layer]) rmin[layer] = rho;
    }

  int layer = ddd->getLayerOffset(); // FIXME: Can be made less stupid
  for (int i = 0; i < disks; ++i) {
    if (countPos[i]) {
      HGCDiskGeomDet* disk = new HGCDiskGeomDet(subdet, +1, layer, zsumPos[i]/countPos[i], rmin[i], rmax[i], radlen_[layer], xi_[layer]);
      addDisk(disk);
    }
    if (countNeg[i]) {
      HGCDiskGeomDet* disk = new HGCDiskGeomDet(subdet, -1, layer, zsumNeg[i]/countNeg[i], rmin[i], rmax[i], radlen_[layer], xi_[layer]);
      addDisk(disk);  
    }
    layer++;
  }
}


template <typename TILES>
void PatternRecognitionbyKF<TILES>::fillHitMap(std::map<DetId,const HGCRecHit*>& hitMap,
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
}

template <typename TILES>
std::vector<TrajectoryMeasurement> PatternRecognitionbyKF<TILES>::measurements(
      const TrajectoryStateOnSurface &tsos, 
      const MeasurementEstimator &mest, 
      const TILES &tiles, 
      int depth){
  
  std::vector<TrajectoryMeasurement> ret;

  // define search window and get bins

  float eta = tsos.globalPosition().eta();
  float phi = tsos.globalPosition().phi();

  float etaMin = eta - etaBinSize;
  float etaMax = eta + etaBinSize;
  float phiMin = phi - phiBinSize;
  float phiMax = phi + phiBinSize;

  auto bins = tiles[depth].searchBoxEtaPhi(etaMin, etaMax, phiMin, phiMax);

  // loop over candidates

  for (int ieta = bins[0]; ieta < bins[1]; ieta++) {
    auto offset = ieta * nPhiBin;
    for (int phi = bins[2]; phi < bins[3]; phi++) {
      int iphi = ((phi % nPhiBin + nPhiBin) % nPhiBin);
      if (!tiles[depth][offset + iphi].empty()) {
        for(auto hit: tiles[depth][offset + iphi]) {
          const auto rec = hitMap.find(hit)->second;
          float energy = rec->energy();

          GlobalPoint globalpoint = rhtools_.getPosition(hit);
          LocalPoint localpoint = tsos.surface().toLocal(globalpoint);

          auto hitptr = std::make_shared<HGCTrackingRecHit>(hit,localpoint,lerr[hit],energy);
          auto mest_pair = mest.estimate(tsos, *hitptr);
          if(mest_pair.first){
            ret.emplace_back(tsos,hitptr,mest_pair.second);
          }
        }
      }
    }
  }
return ret;
}

template <typename TILES>
void PatternRecognitionbyKF<TILES>::makeTracksters(
    const typename PatternRecognitionAlgoBaseT<TILES>::Inputs &input,
    std::vector<Trackster> &result,
    std::unordered_map<int, std::vector<int>> &seedToTracksterAssociation) {}


template <typename TILES>
void PatternRecognitionbyKF<TILES>::init(
    const edm::Event& ev, const edm::EventSetup& es){

    //Get Calo Geometry
    if (es.get<CaloGeometryRecord>().cacheIdentifier() != geomCacheId_) {
      geomCacheId_ = es.get<CaloGeometryRecord>().cacheIdentifier();
      const CaloGeometry* geom = &es.getData(caloGeomToken_);
      rhtools_.setGeometry(*geom);

      makeDisks(8, 26, geom);
      makeDisks(9, 21, geom);
      makeDisks(10,21, geom);

      auto ptrSort = [](const HGCDiskGeomDet *a, const HGCDiskGeomDet *b) -> bool { return (abs(a->position().z())) < (abs(b->position().z())); };
      std::sort(disksPos_.begin(), disksPos_.end(), ptrSort);
      std::sort(disksNeg_.begin(), disksNeg_.end(), ptrSort);
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
    fillHitMap(hitMap, *ee_hits, *fh_hits, *bh_hits);
}


template <typename TILES>
void PatternRecognitionbyKF<TILES>::makeTracksters_verbose(
    const typename PatternRecognitionAlgoBaseT<TILES>::Inputs &input,
    std::vector<KFHit>& kfhits,
    std::vector<KFHit>& prophits,
    float& abs_fail) {

  edm::EventSetup const &es = input.es;
  edm::Event const &ev = input.ev;
  init(ev,es);
  const TILES &tiles = input.tiles;

  // Build tsos from TrackCollection
  edm::Handle<reco::TrackCollection> tracks_h;
  ev.getByToken(trackToken_,tracks_h);

  const reco::TrackCollection& tkx = *tracks_h; 
  if (tkx.empty()){
    std::cout << "No track found!!!! Exited PatternRecognitionbyKF" << std::endl; 
    abs_fail+=1;
    return;
  }

  FreeTrajectoryState fts = trajectoryStateTransform::outerFreeState(tkx.front(),bfield_.product());

  // Propagate through all disks

  // Get first disk

  int zside = fts.momentum().eta() > 0 ? +1 : -1;
  PropagationDirection direction = alongMomentum;
  std::vector<HGCDiskGeomDet*> disks = (zside > 0? disksPos_ : disksNeg_);
  const HGCDiskGeomDet* disk = (zside > 0 ? disksPos_ : disksNeg_).front();
  bool isSilicon = true;

  // Propagation step

  std::vector<TempTrajectory> traj = advanceOneLayer(fts, disk, disks, tiles, direction, isSilicon, TempTrajectory(direction,0));
  if (traj.empty()){
    std::cout << "No track found!!!! Exited PatternRecognitionbyKF" << std::endl; 
    abs_fail+=1;
    return;
  }
  auto lm = traj.back().lastMeasurement();
  TrajectoryStateOnSurface tsos_prop = lm.predictedState();
  KFHit *prophit = new KFHit(tsos_prop, lm.recHit()->geographicalId());
  prophits.push_back(*prophit);

  TrajectoryStateOnSurface tsos_kf = lm.updatedState();
  KFHit *kfhit = new KFHit(tsos_kf, lm.recHit()->geographicalId());
  kfhits.push_back(*kfhit);

  std::vector<TempTrajectory> traj_prop;
  std::vector<TempTrajectory> traj_kf;
  traj_prop.push_back(traj.back());
  traj_kf.push_back(traj.back());

  // Loop over all disks

  unsigned int depth = 2;
  for(disk = nextDisk(disk, direction, disks, isSilicon); disk != nullptr; disk = nextDisk(disk, direction, disks, isSilicon), depth++){
    std::vector<TempTrajectory> newcands_kf;
    for(TempTrajectory & cand : traj_kf){

      TrajectoryStateOnSurface start = cand.lastMeasurement().updatedState();
      std::vector<TempTrajectory> hisTrajs = advanceOneLayer(start, disk, disks, tiles, direction, isSilicon, cand);
      for(TempTrajectory & t : hisTrajs){
        auto lm = t.lastMeasurement();
        newcands_kf.push_back(t);

        TrajectoryStateOnSurface tsos_kf = lm.updatedState();
        KFHit *kfhit = new KFHit(tsos_kf, lm.recHit()->geographicalId());
        kfhits.push_back(*kfhit);
        break;
      }
    }
    traj_kf.swap(newcands_kf);
  }

  // Propagator loop used purely for testing purposes
  direction = alongMomentum;
  disk = (zside > 0 ? disksPos_ : disksNeg_).front();
  isSilicon=1;
  depth = 2;
  for(disk = nextDisk(disk, direction, disks, isSilicon); disk != nullptr; disk = nextDisk(disk, direction, disks, isSilicon), depth++){
    std::vector<TempTrajectory> newcands_prop;
    for(TempTrajectory & cand : traj_prop){

      TrajectoryStateOnSurface start = cand.lastMeasurement().predictedState();
      std::vector<TempTrajectory> hisTrajs = advanceOneLayer(start, disk, disks, tiles, direction, isSilicon, cand);

      for(TempTrajectory & t : hisTrajs){
        auto lm = t.lastMeasurement();
        newcands_prop.push_back(t);
        TrajectoryStateOnSurface tsos_prop = lm.updatedState();
        KFHit *prophit = new KFHit(tsos_prop, lm.recHit()->geographicalId());
        prophits.push_back(*prophit);
        break;

      }
    }
    traj_prop.swap(newcands_prop);
  }
}

template <typename TILES>
void PatternRecognitionbyKF<TILES>::dumpTiles(const TILES &tiles) const {
  std::cout << "Entered dumpTiles" << std::endl;
  constexpr int nEtaBin = TILES::constants_type_t::nEtaBins;
  constexpr int nPhiBin = TILES::constants_type_t::nPhiBins;
  std::cout << nEtaBin << "\t" <<nPhiBin << std::endl;
  auto lastLayerPerSide = static_cast<int>(rhtools_.lastLayer(false));
  std::cout << lastLayerPerSide << std::endl;
  int maxLayer = 2 * lastLayerPerSide - 1;
  int count = 0;
  for (int layer = 0; layer <= maxLayer; layer++) {
    for (int ieta = 0; ieta < nEtaBin; ieta++) {
      auto offset = ieta * nPhiBin;
      for (int phi = 0; phi < nPhiBin; phi++) {
        int iphi = ((phi % nPhiBin + nPhiBin) % nPhiBin);
        if (!tiles[layer][offset + iphi].empty()) {
          std::cout << "Layer: " << layer << " ieta: " << ieta << " phi: " << phi
                                                         << " " << tiles[layer][offset + iphi].size() << std::endl;
          count++;
        }
      }
    }
  }
  std::cout << "Number of RecHits: " << count << std::endl;
}

template <typename TILES>
void PatternRecognitionbyKF<TILES>::fillPSetDescription(edm::ParameterSetDescription &iDesc) {
  iDesc.add<int>("algo_verbosity", 0);
  //iDesc.add<std::string>("propagator", "RungeKuttaTrackerPropagator"); // RungeKutta Propagator
  iDesc.add<std::vector<double>>("radlen",{});
  iDesc.add<std::vector<double>>("xi",{});
  iDesc.add<std::string>("propagator", "PropagatorWithMaterial"); // Analytical Propagator 
  iDesc.add<std::string>("propagatorOpposite", "PropagatorWithMaterialOpposite");
  iDesc.add<std::string>("estimator", "Chi2");
  iDesc.add<edm::InputTag>("tracks", edm::InputTag("generalTracks"));
  iDesc.add<edm::InputTag>("HGCEEInput", edm::InputTag("HGCalRecHit", "HGCEERecHits"));
  iDesc.add<edm::InputTag>("HGCFHInput", edm::InputTag("HGCalRecHit", "HGCHEFRecHits"));
  iDesc.add<edm::InputTag>("HGCBHInput", edm::InputTag("HGCalRecHit", "HGCHEBRecHits"));
}

template class ticl::PatternRecognitionbyKF<TICLLayerTiles>;
template class ticl::PatternRecognitionbyKF<TICLLayerTilesHFNose>;
