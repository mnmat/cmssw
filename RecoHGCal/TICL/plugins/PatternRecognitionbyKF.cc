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

using namespace ticl;

template <typename TILES>
PatternRecognitionbyKF<TILES>::PatternRecognitionbyKF(const edm::ParameterSet &conf, edm::ConsumesCollector iC)
    : PatternRecognitionAlgoBaseT<TILES>(conf, iC),
      caloGeomToken_(iC.esConsumes<CaloGeometry, CaloGeometryRecord>()),
      propName_(conf.getParameter<std::string>("propagator")),
      bfieldtoken_(iC.esConsumes<MagneticField, IdealMagneticFieldRecord>()),
      propagatortoken_(iC.esConsumes<Propagator, TrackingComponentsRecord>(edm::ESInputTag("",propName_))),
      trackToken_(iC.consumes<reco::TrackCollection>(conf.getParameter<edm::InputTag>("tracks"))),
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
const GeomDet * PatternRecognitionbyKF<TILES>::nextDisk(const GeomDet * from, PropagationDirection direction, const std::vector<GeomDet *> &vec) const{
  
  auto it = std::find(vec.begin(), vec.end(), from);
  if (it == vec.end()) throw cms::Exception("LogicError", "nextDisk called with invalid starting disk");
  if (direction == alongMomentum) {
    if(*it == vec.back()) return nullptr;
    return *(++it);    
  } else {
    if (it == vec.begin()) return nullptr;
      return *(--it);
  }
}

template<typename TILES>
void PatternRecognitionbyKF<TILES>::computeAbsorbers(){
  std::map<std::string, float> X0_;
  std::map<std::string, float> lambda_;
  std::map<std::string, float> ZoA_;

  X0_["Fe"] = 13.84;
  X0_["Pb"] = 6.37;
  X0_["Cu"] = 12.86;
  X0_["W"] = 6.76;
  X0_["WCu"] = combineX0(0.75, X0_["W"], 0.25, X0_["Cu"]);
  //X0_["WCu"] = combinedEdX(0.75, X0_["W"], 0.25, X0_["Cu"]);

  lambda_["Fe"] = 132.1;
  lambda_["Pb"] = 199.6;
  lambda_["Cu"] = 137.3;
  lambda_["W"] = 191.9;
  lambda_["WCu"] = combineX0(0.75, lambda_["W"], 0.25, lambda_["Cu"]);
  //lambda_["WCu"] = combinedEdX(0.75, lambda_["W"], 0.25, lambda_["Cu"]);

  ZoA_["Fe"] = 0.466;
  ZoA_["Pb"] = 0.396;
  ZoA_["Cu"] = 0.456;
  ZoA_["W"] = 0.403;
  ZoA_["WCu"] = combinedEdX(0.75, ZoA_["W"], 0.25, ZoA_["Cu"]);

  // see DataFormats/GeometrySurface/interface/MediumProperties.h
  for(auto ij : X0_){
    xi_[ij.first] = X0_[ij.first] * 0.307075 * ZoA_[ij.first] * 0.5;
    std::cout << " " << ij.first << " xi = " << xi_[ij.first] << std::endl;
  }


  //from TDR pag 18
  //radlen = number of X0
  //EE odd layers: 0.748 Pb + 0.068 Fe
  //EE even layers: 0.648 WCu + 0.417 Cu

  //FH (Had first 12) L28: 0.374 Pb + 0.034 Fe + 0.007 Cu + 2.277 Fe + 0.07 Cu
  //FH (Had first 12) L>29: 0.2315 WCu + 1.992 Fe + 0.487 Cu

  //BH (Had last 12): 0.2315 WCu + 3.87 Fe + 0.487 Cu

}


template<typename TILES>
void PatternRecognitionbyKF<TILES>::makeDisks(int subdet, int disks, const CaloGeometry* geom_) {

    // const CaloSubdetectorGeometry *subGeom = subdet < 5 ? geom_->getSubdetectorGeometry(DetId::Forward, subdet) :
    //                                                       geom_->getSubdetectorGeometry(DetId::Hcal, 2);


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

  std::vector<float>  rmax(disks, 0), rmin(disks, 9e9);
  std::vector<double> zsumPos(disks), zsumNeg(disks);
  std::vector<int>    countPos(disks), countNeg(disks);
  const std::vector<DetId> & ids = subGeom->getValidDetIds();
  //if (hgctracking::g_debuglevel > 0) std::cout << "on subdet " << subdet << " I got a total of " << ids.size() << " det ids " << std::endl;

  // This loops over all the valid DetIds of the geometry. This seems like a very roundabout way of getting 47 layers from millions of detids

  for (auto & i : ids) {
    const GlobalPoint & pos = geom_->getPosition(i); 
    float z = pos.z();
    float rho = pos.perp();
    int side = z > 0 ? +1 : -1;

    int layer = std::numeric_limits<unsigned int>::max();
    if (i.det() == DetId::HGCalEE)    layer = HGCSiliconDetId(i).layer() - 1;
    else if(i.det() == DetId::HGCalHSi) layer = HGCSiliconDetId(i).layer() - 1;
    else if (i.det() == DetId::HGCalHSc)  layer = HGCScintillatorDetId(i).layer() - 1;

    (side > 0 ? zsumPos : zsumNeg)[layer] += z;
    (side > 0 ? countPos : countNeg)[layer]++;
    if (rho > rmax[layer]) rmax[layer] = rho;
    if (rho < rmin[layer]) rmin[layer] = rho;
  }
  int j;
  float n = 10.;
  for (int i = 0; i < disks; ++i) {
    float radlen=-1, xi=-1; // see DataFormats/GeometrySurface/interface/MediumProperties.h
    if (materialbudget_ == "Val"){
      if (subdet == 8){ 
        j = 0;
      }
      else if (subdet == 9){
        j = 26;
      }
      radlen = radlen_v[i+j]*n;
      xi = xi_v[i+j]*n;
    }
    else if (materialbudget_ == "AG"){
      switch(subdet) {
        case 8:
          if (i%2 == 0) {
            radlen = 0.748 * xi_["Pb"] + 0.068 * xi_["Fe"] + 0.014 * xi_["Cu"]+n;
            xi = radlen / (0.748 + 0.068 + 0.014) * 1.e-3;
          }
          else{
            radlen = 0.648 * xi_["WCu"] + 0.417 * xi_["Cu"]+n;
            xi = radlen / (0.648 + 0.417) * 1.e-3;
          }
          break;
        case 9:
          if (i == 0){
            radlen = 0.374 * xi_["Pb"] + (0.007+0.07) * xi_["Cu"] + (0.034+2.277) * xi_["Fe"]+n;
            xi = radlen / (0.374 + 0.007+0.07 + 0.034+2.277) * 1.e-3;
          }
          else if(i < 12){
            radlen = 0.2315 * xi_["WCu"] + 0.487 * xi_["Cu"] + 1.992 * xi_["Fe"]+n;
            xi = radlen / (0.2315 + 0.487 + 1.992) * 1.e-3;
          }
          else{
            radlen = 0.2315 * xi_["WCu"] + 0.487 * xi_["Cu"] + 3.870 * xi_["Fe"]+n;
            xi = radlen / (0.2315 + 0.487 + 3.870) * 1.e-3;
          }
          break;
      }
    }

    std::cout << radlen << std::endl;
    std::cout << xi << std::endl;
    if (countPos[i]) {
      //printf("Positive disk %2d at z = %+7.2f   %6.1f <= rho <= %6.1f\n", i+1, zsumPos[i]/countPos[i], rmin[i], rmax[i]);
              //addDisk(new GeomDet(subdet, +1, i+1, zsumPos[i]/countPos[i], rmin[i], rmax[i], radlen, xi));

      GeomDet* disk = new GeomDet(Disk::build(Disk::PositionType(0,0,zsumPos[i]/countPos[i]), Disk::RotationType(), SimpleDiskBounds(rmin[i], rmax[i], -20, 20)).get() );
      if (radlen > 0) {
        (const_cast<Plane &>(disk->surface())).setMediumProperties(MediumProperties(radlen,xi));
      }
      addDisk(disk, 1);
    }
    if (countNeg[i]) {
      //printf("Negative disk %2d at z = %+7.2f   %6.1f <= rho <= %6.1f\n", i+1, zsumNeg[i]/countPos[i], rmin[i], rmax[i]);
      //addDisk(new GeomDet(subdet, -1, i+1, zsumNeg[i]/countNeg[i], rmin[i], rmax[i], radlen, xi));

      GeomDet* disk = new GeomDet(Disk::build(Disk::PositionType(0,0,zsumNeg[i]/countNeg[i]), Disk::RotationType(), SimpleDiskBounds(rmin[i], rmax[i], -20, 20)).get() );
      if (radlen > 0) {
        (const_cast<Plane &>(disk->surface())).setMediumProperties(MediumProperties(radlen,xi));
      }

      addDisk(disk, -1);
    }
  }
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
    std::vector<GlobalPoint> &points,
    std::unordered_map<int, std::vector<int>> &seedToTracksterAssociation) {

  // Initializing (maybe export to own function)

  edm::EventSetup const &es = input.es;
  bfield_ = es.getHandle(bfieldtoken_);
  propagator_ = es.getHandle(propagatortoken_);
  const Propagator &prop = *propagator_;
  const CaloGeometry* geom = &es.getData(caloGeomToken_);

  computeAbsorbers();
  if (disksPos_.size()==0){
    makeDisks(8, 26, geom);
    makeDisks(9, 21, geom);
  }
  std::cout << "Done calcuating Absorbers" <<std::endl;
  // Sort needs to be implemented

  //auto ptrSort = [](const GeomDet *a, const GeomDet *b) -> bool { return (*a) < (*b); };
  //std::sort(disksPos_.begin(), disksPos_.end(), ptrSort);
  //std::sort(disksNeg_.begin(), disksNeg_.end(), ptrSort);

  // Option 1: build from track

  edm::Event const &ev = input.ev;
  edm::Handle<reco::TrackCollection> tracks_h;
  ev.getByToken(trackToken_,tracks_h);
  const reco::TrackCollection& tkx = *tracks_h;  

  FreeTrajectoryState fts = trajectoryStateTransform::outerFreeState(tkx.front(),bfield_.product());
  std::cout << fts.position().x() << std::endl;
  std::cout << fts.position().y() << std::endl;
  std::cout << fts.position().z() << std::endl;

  // Option 2: from SeedingRegion
  // The tracks are chosen as they contain error information which the seedingregion does not.
  // This however means that the propagation to the first layer is done twice: once by the seedingregionproducer and again here in the next step.

  /*
  edm::Event const &ev = input.ev;
  std::cout << "Loaded Event" << std::endl;
  //edm::Ref<reco::TrackCollection> myRef(input.regions.front().collectionID);
  edm::Handle<reco::TrackCollection> tracks_h_seed;
  std::cout<<input.regions.front().collectionID<<std::endl;
  ev.get(input.regions.front().collectionID, tracks_h_seed);
  const reco::TrackCollection& tkx_seed = * tracks_h_seed;
  std::cout << "Loaded TrackCollection Seed" << std::endl;
  
  // Create FTS
  std::cout << "Create FTS" <<std::endl;
  FreeTrajectoryState fts = trajectoryStateTransform::outerFreeState(tkx_seed.front(),bfield_.product());
  std::cout << "Has Error: " << fts.hasError()<<std::endl;
  */

  // Propagate through all disks

  // Get first disk

  int zside = fts.momentum().eta() > 0 ? +1 : -1;
  PropagationDirection direction = alongMomentum;
  std::vector<GeomDet*> disks = (zside > 0? disksPos_ : disksNeg_);
  const GeomDet* disk = (zside > 0 ? disksPos_ : disksNeg_).front();
  //std::cout << "Loaded first disk" <<std::endl;


  // Propagation step

  TrajectoryStateOnSurface tsos = prop.propagate(fts, disk->surface());
  GlobalPoint gp = tsos.globalPosition();
  points.push_back(gp);

  // Loop over all disks

  unsigned int depth = 2;
  for(disk = nextDisk(disk, direction, disks); disk != nullptr; disk = nextDisk(disk, direction, disks), depth++){
    tsos = prop.propagate(tsos, disk->surface());
    points.push_back(tsos.globalPosition());  
  }

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
  iDesc.add<edm::InputTag>("tracks", edm::InputTag("generalTracks"));
  iDesc.add<std::string>("eid_input_name", "input");
  iDesc.add<std::string>("eid_output_name_energy", "output/regressed_energy");
  iDesc.add<std::string>("eid_output_name_id", "output/id_probabilities");
  iDesc.add<double>("eid_min_cluster_energy", 1.);
  iDesc.add<int>("eid_n_layers", 50);
  iDesc.add<int>("eid_n_clusters", 10);
  iDesc.add<std::string>("materialbudget", "Val"); //"Val", "AG", "custom"
}

template class ticl::PatternRecognitionbyKF<TICLLayerTiles>;
template class ticl::PatternRecognitionbyKF<TICLLayerTilesHFNose>;
