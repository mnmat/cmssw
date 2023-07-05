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

    std::vector<HGCDiskGeomDet *> disksSiPos, disksScPos;
		std::vector<HGCDiskGeomDet *> disksSiNeg, disksScNeg;

    makeDisks(hgcalEEId, geom,disksSiPos, disksSiNeg);
    makeDisks(hgcalHSiId, geom,disksSiPos, disksSiNeg);
    makeDisks(hgcalHScId, geom,disksScPos, disksScNeg);

    auto ptrSort = [](const HGCDiskGeomDet *a, const HGCDiskGeomDet *b) -> bool { return (abs(a->position().z())) < (abs(b->position().z())); };
    std::sort(disksSiPos.begin(), disksSiPos.end(), ptrSort);
    std::sort(disksSiNeg.begin(), disksSiNeg.end(), ptrSort);
    std::sort(disksScPos.begin(), disksScPos.end(), ptrSort);
    std::sort(disksScNeg.begin(), disksScNeg.end(), ptrSort);

    offset = disksSiPos.size()-disksScPos.size();
    std::cout << "Exited makeDisks" << std::endl;
    std::cout << disksSiPos.size() << std::endl;
    std::cout << disksSiNeg.size() << std::endl;
    std::cout << disksScPos.size() << std::endl;
    std::cout << disksScNeg.size() << std::endl;

    std::cout << "Sorted makeDisks" << std::endl;
    makeDiskLayers(disksSiPos,disksScPos);
		makeDiskLayers(disksSiNeg,disksScNeg);

    for(auto disk: diskLayerPos_){
      if (disk->first){
        std::cout << disk->first->layer() << ", " << disk->second->layer() << std::endl;
      } else{
        std::cout << disk->second->layer() << std::endl;
      }
    }

}

HGCTracker::~HGCTracker(){}

void HGCTracker::makeDiskLayers(std::vector<HGCDiskGeomDet*>disksSi, std::vector<HGCDiskGeomDet*> disksSc){
  std::cout << "Offset: " << offset  << std::endl;
	for(int layer=0; layer<lastLayer; layer++){
    std::cout << layer << std::endl;
		if (layer<offset){
			addDiskLayer(disksSi[layer]);
		} else{
      std::cout << "Bla" << std::endl;
			addDiskLayer(disksSi[layer],disksSc[layer-offset]);
		}
	}
}

void HGCTracker::makeDisks(int subdet, const CaloGeometry* geom_, std::vector<HGCDiskGeomDet*> &disksPos, std::vector<HGCDiskGeomDet*> &disksNeg) {
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
  if (subdet == DetId::HGCalHSc) offset = ddd->getLayerOffset();
  std::cout << "------------------- Layer Offset: " << layer << ", NumDisks: " << numdisks << std::endl;
  for (int i = 0; i < numdisks; ++i) {
    std::cout << layer << std::endl;
    if (countPos[i]) {
      HGCDiskGeomDet* disk = new HGCDiskGeomDet(subdet, +1, layer+1, zsumPos[i]/countPos[i], rmin[i], rmax[i], radlen_[layer], xi_[layer]);
      disksPos.push_back(disk);
    }
    if (countNeg[i]) {
      HGCDiskGeomDet* disk = new HGCDiskGeomDet(subdet, -1, layer+1, zsumNeg[i]/countNeg[i], rmin[i], rmax[i], radlen_[layer], xi_[layer]);
      disksNeg.push_back(disk);
    }
    layer++;
  }
}

const HGCDiskLayer* HGCTracker::nextDisk(const HGCDiskLayer * from, 
                                          PropagationDirection direction, 
                                          bool isSilicon) const{

  std::cout << "Entered nextDisk" << std::endl;
  //std::cout << from->
  const std::vector<HGCDiskLayer *> & vec = (from->second->zside() > 0 ? diskLayerPos_ : diskLayerNeg_);
  auto it = std::find(vec.begin(), vec.end(),from);
  if (it == vec.end()) throw cms::Exception("LogicError", "nextDisk called with invalid starting disk"); 
  if (direction == alongMomentum){
    if ((*it == vec.back()) || (*it == vec.rbegin()[1])) return nullptr; // if disk last silicon OR last scintillator disk
    std::cout << "exit next Disk" << std::endl;
    return *(++it);
  } else{
    if (it == vec.begin()) return nullptr;
      std::cout << "exit next Disk" << std::endl;
      return *(--it);
  }
  std::cout << "exit next Disk" << std::endl;
  return nullptr; // Returns nullptr if nextDisk reached end of detector
}

#include "FWCore/Utilities/interface/typelookup.h"
TYPELOOKUP_DATA_REG(HGCTracker);
