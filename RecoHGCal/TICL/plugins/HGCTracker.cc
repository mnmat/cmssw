#include "HGCTracker.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"

#include "Geometry/CommonTopologies/interface/HGCDiskGeomDet.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"


HGCTracker::HGCTracker(const CaloGeometry* geom,
                            const std::vector<double> radlen,
                            const std::vector<double> xi):
                            geom_(geom),
                            radlen_(radlen),
                            xi_(xi)  {


    rhtools_.setGeometry(*geom);
    int hgcalEEId = DetId::HGCalEE;
    int hgcalHSiId = DetId::HGCalHSi;
    int hgcalHScId = DetId::HGCalHSc;

    makeDisks(hgcalEEId, geom);
    makeDisks(hgcalHSiId, geom);
    makeDisks(hgcalHScId, geom);

    auto ptrSort = [](const HGCDiskGeomDet *a, const HGCDiskGeomDet *b) -> bool { return (abs(a->position().z())) < (abs(b->position().z())); };
    std::sort(disksPos_.begin(), disksPos_.end(), ptrSort);
    std::sort(disksNeg_.begin(), disksNeg_.end(), ptrSort);
}

HGCTracker::~HGCTracker(){}

void HGCTracker::makeDisks(int subdet, const CaloGeometry* geom_) {
    const CaloSubdetectorGeometry *subGeom = geom_->getSubdetectorGeometry(DetId::Detector(subdet), ForwardSubdetector::ForwardEmpty);
    auto hgcalGeom = static_cast<const HGCalGeometry*>(subGeom);
    const HGCalDDDConstants* ddd = &(hgcalGeom->topology().dddConstants());
    int numdisks = ddd->lastLayer(true);

    std::vector<float>  rmax(numdisks, 0), rmin(numdisks, 9e9);
    std::vector<double> zsumPos(numdisks), zsumNeg(numdisks);
    std::vector<int> countPos(numdisks), countNeg(numdisks);
    const std::vector<DetId> & ids = subGeom->getValidDetIds();
  
    for (auto & i : ids) {
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
  for (int i = 0; i < numdisks; ++i) {
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

const HGCDiskGeomDet* HGCTracker::nextDisk(const HGCDiskGeomDet * from, 
                                          PropagationDirection direction, 
                                          bool isSilicon) const{


  const std::vector<HGCDiskGeomDet *> & vec = (from->zside() > 0 ? disksPos_ : disksNeg_);
  auto it = std::find(vec.begin(), vec.end(),from);
  if (it == vec.end()) throw cms::Exception("LogicError", "nextDisk called with invalid starting disk");
  int currentLayer = (*it)->layer();

  // disks object contains 61 disks corresponding to the z position of the silicon and scintillator layer (has a slight offset).
  // So in the mixed section, a layer has two disks: Silicon and Scintillator. nextDisk() finds the correct disk
  // based on the isSilicon condition by looping through the disks up until the subsequent layer. 
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

const HGCDiskGeomDet* HGCTracker::switchDisk(const HGCDiskGeomDet * from) const{  

  const std::vector<HGCDiskGeomDet *> & vec = (from->zside() > 0 ? disksPos_ : disksNeg_);
                                                         
  auto it = std::find(vec.begin(), vec.end(),from);
  if (it == vec.end()) throw cms::Exception("LogicError", "nextDisk called with invalid starting disk");
  from->isSilicon()? --it: ++it;
  return *(it);
}

#include "FWCore/Utilities/interface/typelookup.h"
TYPELOOKUP_DATA_REG(HGCTracker);
