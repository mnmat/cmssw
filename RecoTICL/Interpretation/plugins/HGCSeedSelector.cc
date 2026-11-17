// HGCSeedSelector.cc
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "DataFormats/Common/interface/RefVector.h"
#include "DataFormats/MuonDetId/interface/GEMDetId.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "Geometry/CommonTopologies/interface/GlobalTrackingGeometry.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

class HGCSeedSelector : public edm::global::EDProducer<> {
public:
  explicit HGCSeedSelector(const edm::ParameterSet& cfg)
      : seedsToken_(consumes(cfg.getParameter<edm::InputTag>("seeds"))),
        geomToken_(esConsumes()),
        bfieldToken_(esConsumes()),
        minPt_(cfg.getParameter<double>("minPt")),
        minAbsEta_(cfg.getParameter<double>("minAbsEta")),
        maxAbsEta_(cfg.getParameter<double>("maxAbsEta")) {
    produces<edm::RefVector<TrajectorySeedCollection>>();
  }

  void produce(edm::StreamID, edm::Event& ev, const edm::EventSetup& es) const override {
    auto out = std::make_unique<edm::RefVector<TrajectorySeedCollection>>();
    const auto seedsH = ev.getHandle(seedsToken_);
    const auto& geom = es.getData(geomToken_);
    const MagneticField* bfield = &es.getData(bfieldToken_);

    for (unsigned int i = 0; i < seedsH->size(); ++i) {
      const TrajectorySeed& seed = (*seedsH)[i];
      const DetId id(seed.startingState().detId());
      if (id.det() != DetId::Muon || id.subdetId() != MuonSubdetId::GEM || GEMDetId(id).station() != 0)
        continue;

      const GeomDet* det = geom.idToDet(id);
      if (det == nullptr)
        continue;
      const auto tsos = trajectoryStateTransform::transientState(seed.startingState(), &det->surface(), bfield);
      if (!tsos.isValid())
        continue;
      const auto p = tsos.globalMomentum();
      const float absEta = std::abs(tsos.globalPosition().eta());
      if (p.perp() < minPt_ || absEta < minAbsEta_ || absEta > maxAbsEta_)
        continue;

      out->push_back(edm::Ref<TrajectorySeedCollection>(seedsH, i));
    }
    ev.put(std::move(out));
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("seeds", edm::InputTag("ancientMuonSeed"));
    desc.add<double>("minPt", 0.);
    desc.add<double>("minAbsEta", 0.);
    desc.add<double>("maxAbsEta", 99.);
    descriptions.add("hgcSeedSelector", desc);
  }

private:
  const edm::EDGetTokenT<TrajectorySeedCollection> seedsToken_;
  const edm::ESGetToken<GlobalTrackingGeometry, GlobalTrackingGeometryRecord> geomToken_;
  const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> bfieldToken_;
  const double minPt_, minAbsEta_, maxAbsEta_;
};

DEFINE_FWK_MODULE(HGCSeedSelector);