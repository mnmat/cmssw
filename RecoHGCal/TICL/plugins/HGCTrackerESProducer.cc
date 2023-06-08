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

  const std::vector<double> radlen_ {0.562719015399016,0.9551903779795146,1.043875789035635,0.8920189121952249,0.9712515382761755,
      0.9638214213266212,1.0131600728257035,0.931610084339171,0.9762915187242845,0.9673823352145134,
      0.9846489838011926,0.9687667855140328,0.9718771132058326,0.9676763827273792,0.9755766684867133,
      0.968786167714319,1.041570725190613,0.905510318469831,1.5448704774827073,0.9671013911948058,
      1.5527918490549046,0.9644685854108038,1.5537767590772382,0.9681776397486216,1.5531742202624652,
      0.9692846975312643,2.8275418030617567,3.003105180494033,3.005570502299292,3.0076357163745446,
      3.005138681398646,3.0106898519208087,3.005449098667002,3.007544959496376,3.0033570538999084,
      3.005449880885637,3.0056943107449103,4.158164083557213,4.1459435055013865,4.087756410389329,
      4.14190322796051,4.143581858605728,4.135287193893763,4.093154468149641,4.1386581126514566,
      4.145881343241816,4.145663423464052};
  const std::vector<double> xi_ {0.0007556128308446044,0.0007601810242278396,0.0007383347392914731,0.0007954300436349306,0.0007079974303217379,
      0.0007587747875111608,0.0007506686955642467,0.000777442404205796,0.0007223198806925819,0.0007536986711079965,
      0.0007687577424509902,0.0007640633376218576,0.0007307246149741174,0.0007518274635026171,0.0007675709426774015,
      0.00075807841348618,0.000743616079227404,0.0007936746547956953,0.0005974640645738568,0.0007508915190685968,
      0.0006382593210233038,0.0007517838017659952,0.0006397287020925379,0.0007579429706626264,0.0006333046151356415,
      0.0007637633509097611,0.0010337707504477352,0.0010108776039545466,0.0010104712217333819,0.0010106368784497104,
      0.0010106003723634636,0.0010158675945647706,0.001010809007077468,0.001009561281210271,0.0010106988945190448,
      0.0010103715066888692,0.0010106495620045056,0.0010210702627516468,0.0010169851234559488,0.0010010856800569403,
      0.0010155575635282914,0.0010146829805806501,0.0010150368006945165,0.0009992234164373278,0.0010145096108228593,
      0.0010126783833862823,0.0010125427722475087};
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
  const CaloGeometry* geom = &(iRecord.get(caloGeomToken_));
  rhtools_.setGeometry(*geom);

  auto product = std::make_unique<HGCTracker>(HGCTracker(geom, radlen_, xi_));
  return product;
}

//define this as a plug-in
DEFINE_FWK_EVENTSETUP_MODULE(HGCTrackerESProducer);