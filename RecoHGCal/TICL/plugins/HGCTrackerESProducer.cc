// -*- C++ -*-
//
// Package:    PlayGround/HGCTrackerESProducer
// Class:      HGCTrackerESProducer
//
/**\class HGCTrackerESProducer

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Mark Nicholas Matthewman
//         Created:  Fri, 19 May 2023 09:58:45 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/Framework/interface/ESProducer.h"

#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"

#include "Geometry/CommonTopologies/interface/HGCDiskGeomDet.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include "HGCTracker.h"


// Need to add #include statements for definitions of
// the data type and record type here

//
// class declaration
//

using namespace edm;

class HGCTrackerESProducer : public edm::ESProducer {
public:
  HGCTrackerESProducer(const edm::ParameterSet&);
  ~HGCTrackerESProducer() override;

  using ReturnType = std::unique_ptr<HGCTracker>;

  ReturnType produce(const CaloGeometryRecord&);

private:

  void makeDisks(int subdet, const CaloGeometry* geom);
  void addDisk(HGCDiskGeomDet *disk) { 
      (disk->zside() > 0 ? disksPos_ : disksNeg_).push_back(disk);
  }
  // ----------member data ---------------------------

  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeomToken_;
  std::string name_ = "Helo, m";

  std::vector<HGCDiskGeomDet *> disksPos_, disksNeg_;
  hgcal::RecHitTools rhtools_;

  const std::vector<double> radlen_ {1.53770787, 0.71064359, 1.45345887, 0.57315113, 1.02882455,
        0.92384098, 1.15461784, 0.72404336, 0.9948446 , 1.0107427 ,
        1.05947235, 0.91730167, 1.20028302, 0.6703572 , 0.98144224,
        1.01024202, 1.08452792, 0.86282149, 1.53923452, 0.99185102,
        1.67874405, 0.70709974, 1.63824099, 0.97162878, 1.74571227,
        0.69011827, 2.92834302, 3.01147101, 3.0583451 , 3.12601533,
        2.85205937, 2.95217992, 3.14263578, 3.07471756, 3.05502943,
        2.82345623, 3.0230636 , 4.29398744, 3.9234094 , 4.27748842,
        3.91229994, 4.23728221, 4.02845205, 4.21537293, 4.32452121,
        3.83363941, 4.32332509};
  const std::vector<double> xi_ {0.00264665, 0.00050171, 0.00081145, 0.0003883 , 0.00049233,
        0.00066116, 0.00059059, 0.00050953, 0.0004874 , 0.00069975,
        0.00051153, 0.00065396, 0.00062511, 0.00046753, 0.00048085,
        0.00070629, 0.00052536, 0.00061608, 0.00068004, 0.0006793 ,
        0.00079585, 0.0005112 , 0.00073527, 0.00067903, 0.00081588,
        0.00049695, 0.00292419, 0.0029901 , 0.00304991, 0.00309913,
        0.00283049, 0.00295079, 0.00310096, 0.00307401, 0.00304154,
        0.00278516, 0.00302603, 0.00425792, 0.00389764, 0.00425274,
        0.0038517 , 0.00426734, 0.00399329, 0.00420363, 0.00429176,
        0.00378711, 0.004295657};
};
//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
HGCTrackerESProducer::HGCTrackerESProducer(const edm::ParameterSet& iConfig) {
  //the following line is needed to tell the framework what
  // data is being produced
  auto cc = setWhatProduced(this);
  caloGeomToken_ = cc.consumes<CaloGeometry>();
  //now do what ever other initialization is needed
}

HGCTrackerESProducer::~HGCTrackerESProducer() {
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

void HGCTrackerESProducer::makeDisks(int subdet, const CaloGeometry* geom_) {

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


// ------------ method called to produce the data  ------------
HGCTrackerESProducer::ReturnType HGCTrackerESProducer::produce(const CaloGeometryRecord& iRecord) {
  std::cout << "Entered HGCTracker produce method" << std::endl;

  const CaloGeometry* geom = &(iRecord.get(caloGeomToken_));
  rhtools_.setGeometry(*geom);

  auto product = std::make_unique<HGCTracker>(HGCTracker(geom, radlen_, xi_));
  std::cout << "Done HGCTracker producer" << std::endl;

  return product;
}

//define this as a plug-in
DEFINE_FWK_EVENTSETUP_MODULE(HGCTrackerESProducer);