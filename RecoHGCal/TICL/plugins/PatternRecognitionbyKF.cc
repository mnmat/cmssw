// Author: Marco Rovere - marco.rovere@cern.ch
// Date: 04/2021
#include <algorithm>
#include <set>
#include <vector>
#include <typeinfo>

#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "PatternRecognitionbyKF.h"

#include "TrackstersPCA.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/TrackReco/interface/Track.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "DataFormats/ForwardDetId/interface/HGCScintillatorDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetId.h"
#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"

#include "DataFormats/GeometrySurface/interface/MediumProperties.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"

#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"

#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimatorBase.h"
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h"
#include "TrackingTools/KalmanUpdators/interface/KFUpdator.h"
#include "TrackingTools/PatternTools/interface/TrajectoryStateUpdator.h"

#include "TrackingTools/PatternTools/interface/TrajMeasLessEstim.h"
#include "DataFormats/TrackingRecHit/interface/InvalidTrackingRecHit.h"

#include "DataFormats/Provenance/interface/Provenance.h"
#include "DataFormats/Provenance/interface/ProcessHistoryRegistry.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"



using namespace ticl;

template <typename TILES>
PatternRecognitionbyKF<TILES>::PatternRecognitionbyKF(const edm::ParameterSet &conf, edm::ConsumesCollector iC)
    : PatternRecognitionAlgoBaseT<TILES>(conf, iC),
      caloGeomToken_(iC.esConsumes<CaloGeometry, CaloGeometryRecord>()),
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
      eidInputName_(conf.getParameter<std::string>("eid_input_name")),
      eidOutputNameEnergy_(conf.getParameter<std::string>("eid_output_name_energy")),
      eidOutputNameId_(conf.getParameter<std::string>("eid_output_name_id")),
      eidMinClusterEnergy_(conf.getParameter<double>("eid_min_cluster_energy")),
      eidNLayers_(conf.getParameter<int>("eid_n_layers")),
      eidNClusters_(conf.getParameter<int>("eid_n_clusters")),
      materialbudget_(conf.getParameter<std::string>("materialbudget")){

      std::cout << propName_ << std::endl;
      std::cout << materialbudget_ <<std::endl;


};

template<typename TILES>
void PatternRecognitionbyKF<TILES>::calculateLocalError(DetId id, const HGCalDDDConstants* ddd){
  if(rhtools_.isSilicon(id)){
    float A;
    if(rhtools_.getSiThickness(id) < 200) A = 1.18; // TODO: replace with non-hardcoded value; hardcoded value from TDR
    else  A = 0.52; // TODO: replace with non-hardcoded value; hardcoded value from TDR
    float a = sqrt(2*A/(3*sqrt(3)));
    double varx = pow(a,4)*5*sqrt(3)/(16*A); // x
    double vary = pow(a,4)*5*sqrt(3)/(16*A); // y 
    lerr[id] = LocalError(varx, 0, vary);
  }
  else{
    const GlobalPoint &pos = rhtools_.getPosition(id);
    double r = sqrt(pos.x()*pos.x() + pos.y()*pos.y());
    auto radiusLayer = ddd->getRadiusLayer(rhtools_.getLayer(id));
    int idx = static_cast<int>(std::lower_bound(radiusLayer.begin(), radiusLayer.end(),r)-radiusLayer.begin());
    float rmax = radiusLayer[idx];
    float rmin = radiusLayer[idx-1];

    double phi = rhtools_.getPhi(id) + M_PI; // radians [0, 2pi]
    double dphi = rhtools_.getScintDEtaDPhi(id).second; // radians
    double phimin = phi - 0.5*dphi;
    double phimax = phi + 0.5*dphi;

    double A = (rmax*rmax - rmin*rmin)*M_PI*dphi/(2*M_PI);

    double ex2 = 1/(8*A) * (pow(rmax,4) - pow(rmin,4)) * (-phimin - sin(phimin)*cos(phimin) + phimax + sin(phimax)*cos(phimax));
    double ex = 1/(3*A) * (pow(rmax,3) - pow(rmin,3)) * (sin(phimax) - sin(phimin));
    double varx = ex2 - ex*ex;

    double ey2 = 1/(8*A) * (pow(rmax,4) - pow(rmin,4)) * (-phimin + sin(phimin)*cos(phimin) + phimax - sin(phimax)*cos(phimax));
    double ey = 1/(3*A) * (pow(rmax,3) - pow(rmin,3)) * (cos(phimin) - cos(phimax));
    double vary = ey2 - ey*ey;

    double varxy = 1/(16*A)*(pow(rmax,4)-pow(rmin,4))*(cos(2*phimin)-cos(2*phimax)) - ex*ey;
    lerr[id] = LocalError(varx, varxy, vary);
    //std::cout << "Phi: " << phi << "\t DPhi: " << dphi << "\t r:" << rmin << "\t R:" << rmax << "\t varx: " << varx << "\t vary: " << vary<< std::endl;

    /*
    std::cout<< ddd->getTrFormN()<<std::endl;
    auto lfb = ddd->getParameter()->layerFrontBH_;
    for(auto it: lfb){
      std::cout << it << std::endl;
    }
    */
    //std::cout<<"First: "<<rhtools_.getScintDEtaDPhi(id).first<<"\t Second: "<<rhtools_.getScintDEtaDPhi(id).second<<std::endl;
    //std::cout<< ddd->getParameter()->scintCellSize(rhtools_.getLayer(id))<<std::endl;
    //std::cout << ddd->getTypeTrap(rhtools_.getLayer(id))<<std::endl;
  }
} 

template<typename TILES>
const HGCDiskGeomDet * PatternRecognitionbyKF<TILES>::switchDisk(const HGCDiskGeomDet * from, const std::vector<HGCDiskGeomDet *> &vec, bool isSilicon) const{
  
  auto it = std::find(vec.begin(), vec.end(),from);
  if (it == vec.end()) throw cms::Exception("LogicError", "nextDisk called with invalid starting disk");
  isSilicon? --it: ++it;
  return *(it);
}


template<typename TILES>
const HGCDiskGeomDet * PatternRecognitionbyKF<TILES>::nextDisk(const HGCDiskGeomDet * from, PropagationDirection direction, const std::vector<HGCDiskGeomDet *> &vec, bool isSilicon) const{

  auto it = std::find(vec.begin(), vec.end(),from);
  if (it == vec.end()) throw cms::Exception("LogicError", "nextDisk called with invalid starting disk");
  int currentLayer = (*it)->layer();
  if (direction == alongMomentum){
    if ((*it == vec.back()) || (*it == vec.rbegin()[1])) return nullptr;
    for(int i = 0; i<3 ; i++){
      if ((*it) == vec.back()) break;
      ++it;
      if(((*(it))->isSilicon() == isSilicon) && ((*(it))->layer()==currentLayer+1)) return *(it);
     }
  } else{
    if (it == vec.begin()) return nullptr;
    for(int i = 0; i<3 ; i++){
      if ((it) == vec.begin()) break;
      --it;
      if(((*(it))->isSilicon() == isSilicon) && ((*(it))->layer()==currentLayer-1)) return *(it);
    }
  }
  std::cout << "Return nullptr" << std::endl;
  return nullptr;
}



template<typename TILES>
template<class Start> 
std::vector<TempTrajectory>
PatternRecognitionbyKF<TILES>::advanceOneLayer(const Start &start, const HGCDiskGeomDet * disk, const std::vector<HGCDiskGeomDet *> &disks, const TILES &tiles, PropagationDirection direction, bool &isSilicon, TempTrajectory traj){

  std::vector<TempTrajectory> ret;
  int depth = disk->layer()+1;
  std::cout << "Layer: " << depth << std::endl;

  const Propagator &prop = (direction == alongMomentum ? *propagator_ : *propagatorOppo_);
  TrajectoryStateOnSurface tsos = prop.propagate(start, disk->surface());
  std::cout << "IsValid: " << tsos.isValid() << std::endl;
  if (!tsos.isValid()) return ret;
  std::cout << "First Propagation done" << std::endl;
  float r = sqrt(pow(tsos.globalPosition().x(),2)+pow(tsos.globalPosition().y(),2));
  if (((disk->rmin() > r) && (!isSilicon)) || (((r > disk->rmax()) && (isSilicon)))) {
    std::cout << "Entered Switchdisk" << std::endl;
    std::cout << "Radius: " << disk->rmin() << "," << r << "," << disk->rmax() << std::endl;
    disk = switchDisk(disk, disks, isSilicon);
    isSilicon = !isSilicon;
    tsos = prop.propagate(start, disk->surface());
    float r = sqrt(pow(tsos.globalPosition().x(),2)+pow(tsos.globalPosition().y(),2));
  }

  // Collect hits with estimate
  depth = disk->layer()+1;
  std::cout << "Layer: " << depth << std::endl;
  auto meas = measurements(tsos, *estimator_, tiles, depth);
  std::sort(meas.begin(), meas.end(),TrajMeasLessEstim());

  if (meas.empty()){
    std::cout << "No measurement found!!!!!!" << std::endl;
  }

  for (const TrajectoryMeasurement &tm : meas){
    TrajectoryStateOnSurface updated = updator_->update(tm.forwardPredictedState(),*tm.recHit());
    ret.push_back(traj.foundHits() ? traj: TempTrajectory(traj.direction(),0));
    ret.back().push(TrajectoryMeasurement(tm.forwardPredictedState(),
                                          updated,
                                          tm.recHit(),
                                          tm.estimate()),
                    tm.estimate());
  }

  if (meas.size() > 0) return ret;

  auto missing = TrackingRecHit::missing;
  ret.push_back(traj.foundHits()? traj : TempTrajectory(traj.direction(),0));
  ret.back().push(TrajectoryMeasurement(tsos, std::make_shared<InvalidTrackingRecHit>(*disk,missing)));

  return ret;
}

template<typename TILES>
void PatternRecognitionbyKF<TILES>::makeDisks(int subdet, int disks, const CaloGeometry* geom_) {

    // const CaloSubdetectorGeometry *subGeom = subdet < 5 ? geom_->getSubdetectorGeometry(DetId::Forward, subdet) :
    //                                                       geom_->getSubdetectorGeometry(DetId::Hcal, 2);



  // FIXME: Load data from xml file or switch between cases

  std::vector<float> radlen_v{0.287578231425533,
    0.652574432663379,
    0.519126219884760,
    0.629329696701607,
    0.470164608944423,
    0.656926010668083,
    0.505405478217971,
    0.645561815316222,
    0.476905911172138,
    0.657129361462176,
    0.494426335907814,
    0.661767816723196,
    0.476697696751714,
    0.656792486637129,
    0.488911438964815,
    0.659428797546396,
    0.519388799886432,
    0.636443052965693,
    0.690440833716593,
    0.656083145894548,
    0.710648981189018,
    0.654738292372854,
    0.711716272998115,
    0.659031423898539,
    0.708941939031308,
    0.661965595922196,
    2.868636061036860,
    2.988866189877523,
    2.989715911253076,
    2.991436303780332,
    2.989599621564430,
    3.002821764306962,
    2.990356629240060,
    2.989573639239810,
    2.988685711688429,
    2.989414475081340,
    2.990114765024607,
    4.164871621966520,
    4.144585482635521,
    4.058009456314032,
    4.137460795188776,
    4.138177409380636,
    4.129951519389949,
    4.060009792683141,
    4.133359022913795,
    4.134560408406477,
    4.13575162757734};
  std::vector<float> xi_v{0.0005110512059408830,
    0.0006831878206769110,
    0.0004973065046027600,
    0.0007055116075429950,
    0.0004840811987581450,
    0.0006815847792258840,
    0.000498840698299919,
    0.0006929527987818480,
    0.0004884872008263560,
    0.0006792860873528980,
    0.0005021346124779450,
    0.0006831033295305010,
    0.0004904917404416280,
    0.0006787315453395390,
    0.0005011512213829390,
    0.0006806752816281460,
    0.0004986591763045000,
    0.0007028556604867630,
    0.0004469247382095310,
    0.0006784016152473840,
    0.0004576588817242630,
    0.0006788591171105650,
    0.00045805568196347000,
    0.0006806926713052890,
    0.00045644714532507900,
    0.0006829423776194960,
    0.0010145335633696400,
    0.0009952585774521000,
    0.00099472493124547,
    0.00099461390470062,
    0.000994829170470436,
    0.0009973866163567710,
    0.000994978298107381,
    0.0009940245879949960,
    0.0009951150189777040,
    0.0009946645572410690,
    0.0009948166566158740,
    0.0010016131009442000,
    0.0009996724454001690,
    0.0009927229142128690,
    0.0009989274416790460,
    0.000998695705935224,
    0.0009987097209326370,
    0.0009919024127419550,
    0.000998719611624487,
    0.0009972693538724180,
    0.0009976091170762653};


    const CaloSubdetectorGeometry *subGeom = geom_->getSubdetectorGeometry(DetId::Detector(subdet), ForwardSubdetector::ForwardEmpty);
    auto geomEE = static_cast<const HGCalGeometry*>(subGeom);
    const HGCalDDDConstants* ddd = &(geomEE->topology().dddConstants());

    std::vector<float>  rmax(disks, 0), rmin(disks, 9e9);
    std::vector<double> zsumPos(disks), zsumNeg(disks);
    std::vector<int> countPos(disks), countNeg(disks);
    const std::vector<DetId> & ids = subGeom->getValidDetIds();
    //if (hgctracking::g_debuglevel > 0) std::cout << "on subdet " << subdet << " I got a total of " << ids.size() << " det ids " << std::endl;
  
    for (auto & i : ids) {

        calculateLocalError(i,ddd);

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


  int layer = ddd->getLayerOffset();
  for (int i = 0; i < disks; ++i) {
    float radlen=-1, xi=-1;
    radlen = radlen_v[layer];
    xi = xi_v[layer];

    if (countPos[i]) {
      //printf("Positive disk %2d at z = %+7.2f   %6.1f <= rho <= %6.1f\n", i+1, zsumPos[i]/countPos[i], rmin[i], rmax[i]);
              //addDisk(new GeomDet(subdet, +1, i+1, zsumPos[i]/countPos[i], rmin[i], rmax[i], radlen, xi));
      HGCDiskGeomDet* disk = new HGCDiskGeomDet(subdet, +1, layer, zsumPos[i]/countPos[i], rmin[i], rmax[i], radlen, xi);
      addDisk(disk);
    }
    if (countNeg[i]) {

      //printf("Negative disk %2d at z = %+7.2f   %6.1f <= rho <= %6.1f\n", i+1, zsumNeg[i]/countPos[i], rmin[i], rmax[i]);
      //addDisk(new GeomDet(subdet, -1, i+1, zsumNeg[i]/countNeg[i], rmin[i], rmax[i], radlen, xi));
      HGCDiskGeomDet* disk = new HGCDiskGeomDet(subdet, -1, layer, zsumNeg[i]/countNeg[i], rmin[i], rmax[i], radlen, xi);
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
} // end of EfficiencyStudies::fillHitMap

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
          //std::cout << "Layer: " << depth << " ieta: " << ieta << " phi: " << phi
          //                                      << "Size: " << tiles[depth][offset + iphi].size() << std::endl;


          auto it = hitMap.find(hit);

          const auto rec = it->second;
          auto detid = rec->detid();

          //typedef TColl::const_iterator const_iterator;
          //auto const_iterator =  it;
          //objref = ref_type(handle, const_iterator)
          GlobalPoint globalpoint = rhtools_.getPosition(it->first);
          std::cout<< "Zpos: " << globalpoint.z() << std::endl;
          LocalPoint localpoint = tsos.surface().toLocal(globalpoint); // 
  
          float energy = rec->energy();


          //std::cout << "Global point: "  << globalpoint << std::endl;
          //std::cout << "Local point: "  << localpoint << "\t" << "Local error:" << localerror << "\t" << "Energy: "<< energy <<  std::endl;


          auto hitptr = std::make_shared<HGCTrackingRecHit>(detid,localpoint,lerr[detid],energy);


        // Test Estimator

          auto mest_pair = mest.estimate(tsos, *hitptr);
          //std::cout << "Est Firs: " << mest_pair.first << "\t" <<"Est Second: " <<  mest_pair.second << std::endl;


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
void PatternRecognitionbyKF<TILES>::makeTracksters_verbose(
    const typename PatternRecognitionAlgoBaseT<TILES>::Inputs &input,
    std::vector<Trackster> &result,
    std::vector<GlobalPoint> &points_kf,
    std::vector<GlobalPoint> &points_prop,
    std::vector<float>& xx_kf,
    std::vector<float>& xy_kf,
    std::vector<float>& yy_kf,
    std::vector<float>& xx_prop,
    std::vector<float>& xy_prop,
    std::vector<float>& yy_prop,
    float& abs_fail,
    std::vector<float>& charge_kf,
    std::unordered_map<int, std::vector<int>> &seedToTracksterAssociation) {

  // Initializing (maybe export to own function)

  edm::EventSetup const &es = input.es;
  edm::Event const &ev = input.ev;

  const TILES &tiles = input.tiles;
  bfield_ = es.getHandle(bfieldtoken_);
  propagator_ = es.getHandle(propagatortoken_);
  propagatorOppo_ = es.getHandle(propagatorOppoToken_);
  estimator_ = es.getHandle(estimatorToken_);
  updator_ = es.getHandle(updatorToken_);
  const CaloGeometry* geom = &es.getData(caloGeomToken_);
  rhtools_.setGeometry(*geom);
  const std::vector<DetId> & ids = geom->getValidDetIds();

  edm::Handle<HGCRecHitCollection> ee_hits;
  edm::Handle<HGCRecHitCollection> fh_hits;
  edm::Handle<HGCRecHitCollection> bh_hits;

  /*

  std::vector<float> zpos;
  std::vector<unsigned int> layers_si;
  std::vector<unsigned int> layers_sci;
  std::vector<unsigned int> layers;
  int count;

  for (auto & i : ids) {
    count++;
    const GlobalPoint & pos = geom->getPosition(i); 
    float z = pos.z();
    if(std::find(zpos.begin(), zpos.end(), z)==zpos.end()){
      zpos.push_back(z);
    }
    unsigned int layer = rhtools_.getLayerWithOffset(i);
    if(std::find(layers.begin(), layers.end(), layer)==layers.end()){
      layers.push_back(layer);
    }


    if(rhtools_.isSilicon(i)){
      if((std::find(layers_si.begin(), layers_si.end(), layer)==layers_si.end()) && (z<0)){
        layers_si.push_back(layer);
      }
    }
    else{
      if(std::find(layers_sci.begin(), layers_sci.end(), layer)==layers_sci.end()){
        layers_sci.push_back(layer);
      }
    }
  }


  std::cout << "Zpositions" << std::endl;

  for(int i=0; i<int(zpos.size());i++){
    std::cout << zpos[i] << std::endl;
  }

  std::cout << "LAyers" << std::endl;

  for(int i=0; i<int(layers.size());i++){
    std::cout << layers[i] << std::endl;
  }

  std::cout << "Number of cells for geometry: " << count << std::endl;
  */

  if (disksPos_.size()==0){
    makeDisks(8, 26, geom);
    makeDisks(9, 21, geom);
    makeDisks(10,21, geom);

    auto ptrSort = [](const HGCDiskGeomDet *a, const HGCDiskGeomDet *b) -> bool { return (abs(a->position().z())) < (abs(b->position().z())); };
    std::sort(disksPos_.begin(), disksPos_.end(), ptrSort);
    std::sort(disksNeg_.begin(), disksNeg_.end(), ptrSort);
  }

 // Option 2: from SeedingRegion
  // The tracks are chosen as they contain error information which the seedingregion does not.
  // This however means that the propagation to the first layer is done twice: once by the seedingregionproducer and again here in the next step.
  /*
  std::cout << "Getting Seed" << std::endl;

  //std::cout<<input.regions.front().collectionID<<std::endl;

  //edm::Ref<reco::TrackCollection> myRef(input.regions.front().collectionID);
  edm::Handle<reco::TrackCollection> tracks_h_seed;
  ev.get(input.regions.front().collectionID, tracks_h_seed);
  const reco::TrackCollection& tkx_seed = * tracks_h_seed;
  
  // Create FTS
  //std::cout << "Create FTS" <<std::endl;
  //FreeTrajectoryState test = trajectoryStateTransform::outerFreeState(tkx_seed.front(),bfield_.product());
  //std::cout << "Has Error: " << test.hasError()<<std::endl;
  //std::cout << "Seeding Region: \t"<< "z: " << test.position().z() << std::endl;
  */
  

  // Option 1: build from track
  edm::Handle<reco::TrackCollection> tracks_h;
  ev.getByToken(trackToken_,tracks_h);
  /*
  const auto& provenance = *tracks_h.provenance();
  const auto& processName = provenance.processName();
  std::cout<<processName<<std::endl;
  std::cout << "branch Name: \t" << provenance.branchName() << std::endl;
  std::cout << "class name: \t" << provenance.className() << std::endl;
  std::cout << "module name: \t" << provenance.moduleName() << std::endl;
  std::cout << "module label: \t" << provenance.moduleLabel() << std::endl;
  std::cout << "product instance name: \t" << provenance.productInstanceName() << std::endl;
  std::cout << "friendly class name: \t" << provenance.friendlyClassName() << std::endl;
  */
  const reco::TrackCollection& tkx = *tracks_h; 
  if (tkx.empty()){
    std::cout << "No track found!!!! Exited PatternRecognitionbyKF" << std::endl; 
    abs_fail+=1;
    return;
  }

  FreeTrajectoryState fts = trajectoryStateTransform::outerFreeState(tkx.front(),bfield_.product());
  std::cout << "x: " << fts.position().x() <<", y:" << fts.position().y()<<", z: " << fts.position().z() << std::endl;

  // Get Rechits

  ev.getByToken(hgcalRecHitsEEToken_, ee_hits);
  ev.getByToken(hgcalRecHitsFHToken_, fh_hits);
  ev.getByToken(hgcalRecHitsBHToken_, bh_hits);
  fillHitMap(hitMap, *ee_hits, *fh_hits, *bh_hits);

  // Propagate through all disks

  // Get first disk

  int zside = fts.momentum().eta() > 0 ? +1 : -1;
  PropagationDirection direction = alongMomentum;
  std::vector<HGCDiskGeomDet*> disks = (zside > 0? disksPos_ : disksNeg_);
  const HGCDiskGeomDet* disk = (zside > 0 ? disksPos_ : disksNeg_).front();
  bool isSilicon = true;

  // Propagation step

  //const Propagator &prop = (direction == alongMomentum ? *propagator_ : *propagatorOppo_); // FIXME; inconvenient to have to define the propagator here and in the advanceOneLAyer
  //TrajectoryStateOnSurface tsos = prop.propagate(fts, disk->surface());
  //GlobalPoint gp = tsos.globalPosition();

  std::cout << "Advance One Layer" << std::endl;
  std::vector<TempTrajectory> traj = advanceOneLayer(fts, disk, disks, tiles, direction, isSilicon, TempTrajectory(direction,0));
  if (traj.empty()){
    std::cout << "No track found!!!! Exited PatternRecognitionbyKF" << std::endl; 
    abs_fail+=1;
    return;
  }
  std::cout << "Exited Advance One Layer" << std::endl;
  TrajectoryStateOnSurface tsos_prop = traj.back().lastMeasurement().predictedState();
  points_prop.push_back(tsos_prop.globalPosition());
  xx_prop.push_back(tsos_prop.localError().positionError().xx());
  xy_prop.push_back(tsos_prop.localError().positionError().xy());
  yy_prop.push_back(tsos_prop.localError().positionError().yy());
  TrajectoryStateOnSurface tsos_kf = traj.back().lastMeasurement().updatedState();
  points_kf.push_back(tsos_kf.globalPosition());
  xx_kf.push_back(tsos_kf.localError().positionError().xx());
  xy_kf.push_back(tsos_kf.localError().positionError().xy());
  yy_kf.push_back(tsos_kf.localError().positionError().yy());


  std::vector<TempTrajectory> traj_prop;
  std::vector<TempTrajectory> traj_kf;
  traj_prop.push_back(traj.back());
  traj_kf.push_back(traj.back());

  // Loop over all disks

  std::cout << "Start KF" << std::endl;

  unsigned int depth = 2;
  for(disk = nextDisk(disk, direction, disks, isSilicon); disk != nullptr; disk = nextDisk(disk, direction, disks, isSilicon), depth++){
    std::vector<TempTrajectory> newcands_kf;
    for(TempTrajectory & cand : traj_kf){

      TrajectoryStateOnSurface start = cand.lastMeasurement().updatedState();
      std::vector<TempTrajectory> hisTrajs = advanceOneLayer(start, disk, disks, tiles, direction, isSilicon, cand);
      for(TempTrajectory & t : hisTrajs){
        charge_kf.push_back(t.lastMeasurement().updatedState().charge());
        newcands_kf.push_back(t);
        points_kf.push_back(t.lastMeasurement().updatedState().globalPosition());
        xx_kf.push_back(t.lastMeasurement().updatedState().localError().positionError().xx());
        xy_kf.push_back(t.lastMeasurement().updatedState().localError().positionError().xy());
        yy_kf.push_back(t.lastMeasurement().updatedState().localError().positionError().yy());
        break;
      }
    }
    traj_kf.swap(newcands_kf);
  }

  /*

  // Reverse direction

  direction = oppositeToMomentum;
  disk = (isSilicon ? disks.end()[-1] : disks.end()[-2]);


  for(disk = nextDisk(disk, direction, disks, isSilicon); disk != nullptr; disk = nextDisk(disk, direction, disks, isSilicon), depth++){
    if(disk->layer()==33) isSilicon=1;
    std::vector<TempTrajectory> newcands_kf;
    for(TempTrajectory & cand : traj_kf){

      TrajectoryStateOnSurface start = cand.lastMeasurement().updatedState();
      std::vector<TempTrajectory> hisTrajs = advanceOneLayer(start, disk, disks, tiles, direction, isSilicon, cand);

      for(TempTrajectory & t : hisTrajs){
        newcands_kf.push_back(t);
        points_kf.push_back(t.lastMeasurement().updatedState().globalPosition());    
        xx_kf.push_back(t.lastMeasurement().updatedState().localError().positionError().xx());
        xy_kf.push_back(t.lastMeasurement().updatedState().localError().positionError().xy());
        yy_kf.push_back(t.lastMeasurement().updatedState().localError().positionError().yy());


        break;
      }
    }
    traj_kf.swap(newcands_kf);
  }

  */

  std::cout << "Start Propagator" << std::endl;

  direction = alongMomentum;
  disk = (zside > 0 ? disksPos_ : disksNeg_).front();
  isSilicon=1;
  depth = 2;
  for(disk = nextDisk(disk, direction, disks, isSilicon); disk != nullptr; disk = nextDisk(disk, direction, disks, isSilicon), depth++){
    std::vector<TempTrajectory> newcands_prop;
    for(TempTrajectory & cand : traj_prop){

      TrajectoryStateOnSurface start = cand.lastMeasurement().predictedState();
      std::cout <<"Propagator:" << start.globalPosition().z() << std::endl;
      std::vector<TempTrajectory> hisTrajs = advanceOneLayer(start, disk, disks, tiles, direction, isSilicon, cand);

      for(TempTrajectory & t : hisTrajs){
        newcands_prop.push_back(t);
        points_prop.push_back(t.lastMeasurement().predictedState().globalPosition());        
        xx_prop.push_back(t.lastMeasurement().predictedState().localError().positionError().xx());
        xy_prop.push_back(t.lastMeasurement().predictedState().localError().positionError().yy());
        yy_prop.push_back(t.lastMeasurement().predictedState().localError().positionError().xy());
        
        break;
      }
    }
    traj_prop.swap(newcands_prop);
  }

  /*

  // Reverse direction

  direction = oppositeToMomentum;
  disk = (isSilicon ? disks.end()[-1] : disks.end()[-2]);


  for(disk = nextDisk(disk, direction, disks, isSilicon); disk != nullptr; disk = nextDisk(disk, direction, disks, isSilicon), depth++){
    if(disk->layer()==33) isSilicon=1;

    std::vector<TempTrajectory> newcands_prop;

    for(TempTrajectory & cand : traj_prop){

      TrajectoryStateOnSurface start = cand.lastMeasurement().predictedState();
      std::vector<TempTrajectory> hisTrajs = advanceOneLayer(start, disk, disks, tiles, direction, isSilicon, cand);

      for(TempTrajectory & t : hisTrajs){
        newcands_prop.push_back(t);
        points_prop.push_back(t.lastMeasurement().predictedState().globalPosition());        
        xx_prop.push_back(t.lastMeasurement().predictedState().localError().positionError().xx());
        xy_prop.push_back(t.lastMeasurement().predictedState().localError().positionError().yy());
        yy_prop.push_back(t.lastMeasurement().predictedState().localError().positionError().xy());
        break;
      }
    }
    traj_prop.swap(newcands_prop);
  }
  */

  /*
  for(int i = 0; i != 47; i++ ){
    std::cout << points[i].x()<<","<<points[i].y()<<","<<points[i].z() << "," << points[points.size()-i-1].x()<< "," << points[points.size()-i-1].y() <<"," << points[points.size()-i-1].z() <<std::endl;
  }
  */
  /*

  if (input.regions.empty())
    return;

  const int eventNumber = input.ev.eventAuxiliary().event();
  if (PatternRecognitionAlgoBaseT<TILES>::algo_verbosity_ > PatternRecognitionAlgoBaseT<TILES>::Advanced) {
    edm::LogVerbatim("PatternRecognitionbyKFD") << "New Event";
  }

  edm::EventSetup const &es = input.es;
  const CaloGeometry &geom = es.getData(caloGeomToken_);
  rhtools_.setGeometry(geom);

  clusters_.clear();
  clusters_.resize(2 * rhtools_.lastLayer(false));
  std::vector<std::pair<int, int>> layerIdx2layerandSoa;  //used everywhere also to propagate cluster masking

  layerIdx2layerandSoa.reserve(input.layerClusters.size());
  unsigned int layerIdx = 0;
  for (auto const &lc : input.layerClusters) {
    if (input.mask[layerIdx] == 0.) {
      if (PatternRecognitionAlgoBaseT<TILES>::algo_verbosity_ > PatternRecognitionAlgoBaseT<TILES>::Advanced) {
        edm::LogVerbatim("PatternRecognitionbyKFD") << "Skipping masked clustrer: " << layerIdx;
      }
      layerIdx2layerandSoa.emplace_back(-1, -1);
      layerIdx++;
      continue;
    }
    const auto firstHitDetId = lc.hitsAndFractions()[0].first;
    int layer = rhtools_.getLayerWithOffset(firstHitDetId) - 1 +
                rhtools_.lastLayer(false) * ((rhtools_.zside(firstHitDetId) + 1) >> 1);
    assert(layer >= 0);

    layerIdx2layerandSoa.emplace_back(layer, clusters_[layer].x.size());
    float sum_x = 0.;
    float sum_y = 0.;
    float sum_sqr_x = 0.;
    float sum_sqr_y = 0.;
    float ref_x = lc.x();
    float ref_y = lc.y();
    float invClsize = 1. / lc.hitsAndFractions().size();
    for (auto const &hitsAndFractions : lc.hitsAndFractions()) {
      auto const &point = rhtools_.getPosition(hitsAndFractions.first);
      sum_x += point.x() - ref_x;
      sum_sqr_x += (point.x() - ref_x) * (point.x() - ref_x);
      sum_y += point.y() - ref_y;
      sum_sqr_y += (point.y() - ref_y) * (point.y() - ref_y);
    }
    // The variance of X for X uniform in circle of radius R is R^2/4,
    // therefore we multiply the sqrt(var) by 2 to have a rough estimate of the
    // radius. On the other hand, while averaging the x and y radius, we would
    // end up dividing by 2. Hence we omit the value here and in the average
    // below, too.
    float radius_x = sqrt((sum_sqr_x - (sum_x * sum_x) * invClsize) * invClsize);
    float radius_y = sqrt((sum_sqr_y - (sum_y * sum_y) * invClsize) * invClsize);
    clusters_[layer].x.emplace_back(lc.x());
    clusters_[layer].y.emplace_back(lc.y());
    clusters_[layer].radius.emplace_back(radius_x + radius_y);
    clusters_[layer].eta.emplace_back(lc.eta());
    clusters_[layer].phi.emplace_back(lc.phi());
    clusters_[layer].cells.push_back(lc.hitsAndFractions().size());
    clusters_[layer].energy.emplace_back(lc.energy());
    clusters_[layer].isSeed.push_back(false);
    clusters_[layer].clusterIndex.emplace_back(-1);
    clusters_[layer].layerClusterOriginalIdx.emplace_back(layerIdx++);
    clusters_[layer].nearestHigher.emplace_back(-1, -1);
    clusters_[layer].rho.emplace_back(0.f);
    clusters_[layer].delta.emplace_back(std::numeric_limits<float>::max());
  }
  for (unsigned int layer = 0; layer < clusters_.size(); layer++) {
    clusters_[layer].followers.resize(clusters_[layer].x.size());
  }

  auto lastLayerPerSide = static_cast<int>(rhtools_.lastLayer(false));
  int maxLayer = 2 * lastLayerPerSide - 1;
  std::vector<int> numberOfClustersPerLayer(maxLayer, 0);
  for (int i = 0; i <= maxLayer; i++) {
    calculateLocalDensity(input.tiles, i, layerIdx2layerandSoa);
  }
  for (int i = 0; i <= maxLayer; i++) {
    calculateDistanceToHigher(input.tiles, i, layerIdx2layerandSoa);
  }

  auto nTracksters = findAndAssignTracksters(input.tiles, layerIdx2layerandSoa);
  if (PatternRecognitionAlgoBaseT<TILES>::algo_verbosity_ > PatternRecognitionAlgoBaseT<TILES>::Advanced) {
    edm::LogVerbatim("PatternRecognitionbyKFD") << "Reconstructed " << nTracksters << " tracksters" << std::endl;
    dumpClusters(layerIdx2layerandSoa, eventNumber);
  }

  // Build Trackster
  result.resize(nTracksters);

  for (unsigned int layer = 0; layer < clusters_.size(); ++layer) {
    const auto &thisLayer = clusters_[layer];
    if (PatternRecognitionAlgoBaseT<TILES>::algo_verbosity_ > PatternRecognitionAlgoBaseT<TILES>::Advanced) {
      edm::LogVerbatim("PatternRecognitionbyKFD") << "Examining Layer: " << layer;
    }
    for (unsigned int lc = 0; lc < thisLayer.x.size(); ++lc) {
      if (PatternRecognitionAlgoBaseT<TILES>::algo_verbosity_ > PatternRecognitionAlgoBaseT<TILES>::Advanced) {
        edm::LogVerbatim("PatternRecognitionbyKFD") << "Trackster " << thisLayer.clusterIndex[lc];
      }
      if (thisLayer.clusterIndex[lc] >= 0) {
        if (PatternRecognitionAlgoBaseT<TILES>::algo_verbosity_ > PatternRecognitionAlgoBaseT<TILES>::Advanced) {
          edm::LogVerbatim("PatternRecognitionbyKFD") << " adding lcIdx: " << thisLayer.layerClusterOriginalIdx[lc];
        }
        //result[layer].vertices().push_back(1);
        //result[layer].vertex_multiplicity().push_back(1);
        // loop over followers
        
        for (auto [follower_lyrIdx, follower_soaIdx] : thisLayer.followers[lc]) {
          std::array<unsigned int, 2> edge = {
              {(unsigned int)thisLayer.layerClusterOriginalIdx[lc],
               (unsigned int)clusters_[follower_lyrIdx].layerClusterOriginalIdx[follower_soaIdx]}};
          result[thisLayer.clusterIndex[lc]].edges().push_back(edge);

        }
      }
    }
  }

  result.erase(
      std::remove_if(std::begin(result),
                     std::end(result),
                     [&](auto const &v) { return static_cast<int>(v.vertices().size()) < minNumLayerCluster_; }),
      result.end());
  result.shrink_to_fit();

  ticl::assignPCAtoTracksters(result,
                              input.layerClusters,
                              input.layerClustersTime,
                              rhtools_.getPositionLayer(rhtools_.lastLayerEE(false), false).z());

  // run energy regression and ID
  energyRegressionAndID(input.layerClusters, input.tfSession, result);
  if (PatternRecognitionAlgoBaseT<TILES>::algo_verbosity_ > PatternRecognitionAlgoBaseT<TILES>::Advanced) {
    for (auto const &t : result) {
      edm::LogVerbatim("PatternRecognitionbyKFD") << "Barycenter: " << t.barycenter();
      edm::LogVerbatim("PatternRecognitionbyKFD") << "LCs: " << t.vertices().size();
      edm::LogVerbatim("PatternRecognitionbyKFD") << "Energy: " << t.raw_energy();
      edm::LogVerbatim("PatternRecognitionbyKFD") << "Regressed: " << t.regressed_energy();
    }
  }

  // Dump Tracksters information
  if (PatternRecognitionAlgoBaseT<TILES>::algo_verbosity_ > PatternRecognitionAlgoBaseT<TILES>::Advanced) {
    dumpTracksters(layerIdx2layerandSoa, eventNumber, result);
  }

  // Reset internal clusters_ structure of array for next event
  reset();
  if (PatternRecognitionAlgoBaseT<TILES>::algo_verbosity_ > PatternRecognitionAlgoBaseT<TILES>::Advanced) {
    edm::LogVerbatim("PatternRecognitionbyKFD") << std::endl;
  }

  */

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
void PatternRecognitionbyKF<TILES>::energyRegressionAndID(const std::vector<reco::CaloCluster> &layerClusters,
                                                              const tensorflow::Session *eidSession,
                                                              std::vector<Trackster> &tracksters) {
  // Energy regression and particle identification strategy:
  //
  // 1. Set default values for regressed energy and particle id for each trackster.
  // 2. Store indices of tracksters whose total sum of cluster energies is above the
  //    eidMinClusterEnergy_ (GeV) treshold. Inference is not applied for soft tracksters.
  // 3. When no trackster passes the selection, return.
  // 4. Create input and output tensors. The batch dimension is determined by the number of
  //    selected tracksters.
  // 5. Fill input tensors with layer cluster features. Per layer, clusters are ordered descending
  //    by energy. Given that tensor data is contiguous in memory, we can use pointer arithmetic to
  //    fill values, even with batching.
  // 6. Zero-fill features for empty clusters in each layer.
  // 7. Batched inference.
  // 8. Assign the regressed energy and id probabilities to each trackster.
  //
  // Indices used throughout this method:
  // i -> batch element / trackster
  // j -> layer
  // k -> cluster
  // l -> feature

  // set default values per trackster, determine if the cluster energy threshold is passed,
  // and store indices of hard tracksters
  std::vector<int> tracksterIndices;
  std::cout<<"energyRegressionAndID" <<std::endl;
  for (int i = 0; i < static_cast<int>(tracksters.size()); i++) {
    // calculate the cluster energy sum (2)
    // note: after the loop, sumClusterEnergy might be just above the threshold which is enough to
    // decide whether to run inference for the trackster or not
    float sumClusterEnergy = 0.;
    std::cout << "Regressed Energy" << tracksters[i].regressed_energy() <<std::endl;
    for (const unsigned int &vertex : tracksters[i].vertices()) {
      sumClusterEnergy += static_cast<float>(layerClusters[vertex].energy());
      // there might be many clusters, so try to stop early
      if (sumClusterEnergy >= eidMinClusterEnergy_) {
        // set default values (1)
        tracksters[i].setRegressedEnergy(0.f);
        tracksters[i].zeroProbabilities();
        tracksterIndices.push_back(i);
        break;
      }
    }
  }

  std::cout << "Regressed Energy" << tracksters.size() <<std::endl;

  // do nothing when no trackster passes the selection (3)
  int batchSize = static_cast<int>(tracksterIndices.size());
  if (batchSize == 0) {
    return;
  }

  // create input and output tensors (4)
  tensorflow::TensorShape shape({batchSize, eidNLayers_, eidNClusters_, eidNFeatures_});
  tensorflow::Tensor input(tensorflow::DT_FLOAT, shape);
  tensorflow::NamedTensorList inputList = {{eidInputName_, input}};

  std::vector<tensorflow::Tensor> outputs;
  std::vector<std::string> outputNames;
  if (!eidOutputNameEnergy_.empty()) {
    outputNames.push_back(eidOutputNameEnergy_);
  }
  if (!eidOutputNameId_.empty()) {
    outputNames.push_back(eidOutputNameId_);
  }

  // fill input tensor (5)
  for (int i = 0; i < batchSize; i++) {
    const Trackster &trackster = tracksters[tracksterIndices[i]];

    // per layer, we only consider the first eidNClusters_ clusters in terms of energy, so in order
    // to avoid creating large / nested structures to do the sorting for an unknown number of total
    // clusters, create a sorted list of layer cluster indices to keep track of the filled clusters
    std::vector<int> clusterIndices(trackster.vertices().size());
    for (int k = 0; k < (int)trackster.vertices().size(); k++) {
      clusterIndices[k] = k;
    }
    sort(clusterIndices.begin(), clusterIndices.end(), [&layerClusters, &trackster](const int &a, const int &b) {
      return layerClusters[trackster.vertices(a)].energy() > layerClusters[trackster.vertices(b)].energy();
    });

    // keep track of the number of seen clusters per layer
    std::vector<int> seenClusters(eidNLayers_);

    // loop through clusters by descending energy
    for (const int &k : clusterIndices) {
      // get features per layer and cluster and store the values directly in the input tensor
      const reco::CaloCluster &cluster = layerClusters[trackster.vertices(k)];
      int j = rhtools_.getLayerWithOffset(cluster.hitsAndFractions()[0].first) - 1;
      if (j < eidNLayers_ && seenClusters[j] < eidNClusters_) {
        // get the pointer to the first feature value for the current batch, layer and cluster
        float *features = &input.tensor<float, 4>()(i, j, seenClusters[j], 0);

        // fill features
        *(features++) = float(cluster.energy() / float(trackster.vertex_multiplicity(k)));
        *(features++) = float(std::abs(cluster.eta()));
        *(features) = float(cluster.phi());

        // increment seen clusters
        seenClusters[j]++;
      }
    }

    // zero-fill features of empty clusters in each layer (6)
    for (int j = 0; j < eidNLayers_; j++) {
      for (int k = seenClusters[j]; k < eidNClusters_; k++) {
        float *features = &input.tensor<float, 4>()(i, j, k, 0);
        for (int l = 0; l < eidNFeatures_; l++) {
          *(features++) = 0.f;
        }
      }
    }
  }

  // run the inference (7)
  tensorflow::run(const_cast<tensorflow::Session *>(eidSession), inputList, outputNames, &outputs);

  // store regressed energy per trackster (8)

// get the pointer to the energy tensor, dimension is batch x 1
	//float *energy = outputs[0].flat<float>().data();


  float *energy;
  float e = 654321.0f;
  energy = &e;
  for (const int &i : tracksterIndices) {
    tracksters[i].setRegressedEnergy(*(energy++));
  }

  // store id probabilities per trackster (8)
  // get the pointer to the id probability tensor, dimension is batch x id_probabilities.size()
  int probsIdx = eidOutputNameEnergy_.empty() ? 0 : 1;
  //float *probs = outputs[probsIdx].flat<float>().data();
  float *probs;
  float val=0.42f;
  probs=&val;

  for (const int &i : tracksterIndices) {
    tracksters[i].setProbabilities(probs);
    probs += tracksters[i].id_probabilities().size();
  }
  std::cout<<"End of PID and EREG" <<std::endl;

}

template <typename TILES>
void PatternRecognitionbyKF<TILES>::fillPSetDescription(edm::ParameterSetDescription &iDesc) {
  iDesc.add<int>("algo_verbosity", 0);
  //iDesc.add<std::string>("propagator", "RungeKuttaTrackerPropagator"); //PropagatorWithMaterial
  iDesc.add<std::string>("propagator", "PropagatorWithMaterial"); //PropagatorWithMaterial
  iDesc.add<std::string>("propagatorOpposite", "PropagatorWithMaterialOpposite");
  iDesc.add<std::string>("estimator", "Chi2");
  iDesc.add<edm::InputTag>("tracks", edm::InputTag("generalTracks"));
  iDesc.add<std::string>("eid_input_name", "input");
  iDesc.add<std::string>("eid_output_name_energy", "output/regressed_energy");
  iDesc.add<std::string>("eid_output_name_id", "output/id_probabilities");
  iDesc.add<double>("eid_min_cluster_energy", 1.);
  iDesc.add<int>("eid_n_layers", 50);
  iDesc.add<int>("eid_n_clusters", 10);
  iDesc.add<std::string>("materialbudget", "Val"); //"Val", "AG", "custom"
  iDesc.add<edm::InputTag>("HGCEEInput", edm::InputTag("HGCalRecHit", "HGCEERecHits"));
  iDesc.add<edm::InputTag>("HGCFHInput", edm::InputTag("HGCalRecHit", "HGCHEFRecHits"));
  iDesc.add<edm::InputTag>("HGCBHInput", edm::InputTag("HGCalRecHit", "HGCHEBRecHits"));
}

template class ticl::PatternRecognitionbyKF<TICLLayerTiles>;
template class ticl::PatternRecognitionbyKF<TICLLayerTilesHFNose>;
