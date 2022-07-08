#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "DataFormats/ForwardDetId/interface/HGCScintillatorDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetId.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"


#include <vector>

// Analyzer to test the LocalError (Variance and Covariance) calculations found in the recHitTools for various cell types

class HGCalToolTesterLocalError : public edm::one::EDAnalyzer<> {
public:
  explicit HGCalToolTesterLocalError(const edm::ParameterSet&);
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  void beginJob() override {}
  void analyze(edm::Event const& iEvent, edm::EventSetup const&) override;
  void endJob() override {}

private:
  const edm::ESGetToken<CaloGeometry, CaloGeometryRecord> tok_geom_;
  hgcal::RecHitTools tool_;
  std::vector<DetId> si120, si200, si300, sc;
};

HGCalToolTesterLocalError::HGCalToolTesterLocalError(const edm::ParameterSet&)
    : tok_geom_(esConsumes<CaloGeometry, CaloGeometryRecord>()) {}

void HGCalToolTesterLocalError::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  descriptions.add("hgcalToolTesterLocalError", desc);
}

void HGCalToolTesterLocalError::analyze(const edm::Event& /*iEvent*/, const edm::EventSetup& iSetup) {
  const CaloGeometry geo = iSetup.getData(tok_geom_);
  tool_.setGeometry(geo);
  std::vector<DetId::Detector> dets = {DetId::HGCalEE, DetId::HGCalHSi, DetId::HGCalHSc};
  for (const auto& subdet : dets) {
    const CaloSubdetectorGeometry *subGeom = geo.getSubdetectorGeometry(DetId::Detector(subdet), ForwardSubdetector::ForwardEmpty);
    const std::vector<DetId> & ids = subGeom->getValidDetIds();
    // Todo: Find better way of getting the different detids
    for (const auto& id : ids){
      if(tool_.isSilicon(id)){
        if(119. < tool_.getSiThickness(id) && tool_.getSiThickness(id) < 121.){
          si120.push_back(id);
        }
        else if(199. < tool_.getSiThickness(id) && tool_.getSiThickness(id) < 201.){
          si200.push_back(id);
        }
        else if(299. < tool_.getSiThickness(id) && tool_.getSiThickness(id) < 301.){
          si300.push_back(id);
        }
      }
      else{
        sc.push_back(id);
      }
    }
  }
  edm::LogVerbatim("HGCalGeom") << "Local Error of Silicon cells with 120µm thickness: xx=" << tool_.getLocalError(si120[0]).xx()
                                << ", xy=" << tool_.getLocalError(si120[0]).xy() << ", yy=" << tool_.getLocalError(si120[0]).yy();
  edm::LogVerbatim("HGCalGeom") << "Local Error of Silicon cells with 200µm thickness: xx=" << tool_.getLocalError(si200[0]).xx()
                                << ", xy=" << tool_.getLocalError(si200[0]).xy() << ", yy=" << tool_.getLocalError(si200[0]).yy();   
  edm::LogVerbatim("HGCalGeom") << "Local Error of Silicon cells with 300µm thicknes: xx=" << tool_.getLocalError(si300[0]).xx()
                                << ", xy=" << tool_.getLocalError(si300[0]).xy() << ", yy=" << tool_.getLocalError(si300[0]).yy(); 
  edm::LogVerbatim("HGCalGeom") << "Local Error of Scintillator cell with DetID " << sc[0].rawId() << ", with an eta-phi span [" << tool_.getScintDEtaDPhi(sc[0]).first 
                                << ", " << tool_.getScintDEtaDPhi(sc[0]).second << "]: xx=" << tool_.getLocalError(sc[0]).xx() 
                                << ", xy=" << tool_.getLocalError(sc[0]).xy() << ", yy=" << tool_.getLocalError(sc[0]).yy();
}

DEFINE_FWK_MODULE(HGCalToolTesterLocalError);
