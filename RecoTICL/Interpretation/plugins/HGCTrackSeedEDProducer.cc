#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "Geometry/CommonTopologies/interface/GlobalTrackingGeometry.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

class HGCTrackSeedEDProducer : public edm::global::EDProducer<> {
public:
  explicit HGCTrackSeedEDProducer(const edm::ParameterSet&);
  static void fillDescriptions(edm::ConfigurationDescriptions&);
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

private:
  const edm::EDGetTokenT<reco::TrackCollection> tracksToken_;
  const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> bfieldToken_;
  const edm::ESGetToken<GlobalTrackingGeometry, GlobalTrackingGeometryRecord> geomToken_;
  const StringCutObjectSelector<reco::Track> cut_;
  const PropagationDirection direction_;
};

HGCTrackSeedEDProducer::HGCTrackSeedEDProducer(const edm::ParameterSet& cfg)
    : tracksToken_(consumes(cfg.getParameter<edm::InputTag>("InputCollection"))),
      bfieldToken_(esConsumes()),
      geomToken_(esConsumes()),
      cut_(cfg.getParameter<std::string>("cut")),
      direction_(cfg.getParameter<bool>("alongMomentum") ? alongMomentum : oppositeToMomentum) {
  produces<TrajectorySeedCollection>();
  produces<reco::TrackRefVector>("sourceTracks");  // sourceTracks[i] <-> seeds[i]
}

void HGCTrackSeedEDProducer::produce(edm::StreamID, edm::Event& ev, const edm::EventSetup& es) const {
  auto seeds = std::make_unique<TrajectorySeedCollection>();
  auto sources = std::make_unique<reco::TrackRefVector>();

  const auto tracks = ev.getHandle(tracksToken_);
  const MagneticField* bfield = &es.getData(bfieldToken_);
  const GlobalTrackingGeometry& geom = es.getData(geomToken_);

  for (unsigned int i = 0; i < tracks->size(); ++i) {
    const reco::Track& tk = (*tracks)[i];
    if (!cut_(tk))
      continue;

    TrajectoryStateOnSurface tsos = trajectoryStateTransform::outerStateOnSurface(tk, geom, bfield);
    if (!tsos.isValid()) {
      edm::LogWarning("HGCTrackSeedEDProducer") << "Invalid outer state for track " << i << ", skipped.";
      continue;
    }

    seeds->emplace_back(trajectoryStateTransform::persistentState(tsos, tk.outerDetId()),
                        TrajectorySeed::RecHitContainer{},
                        direction_);
    sources->push_back(reco::TrackRef(tracks, i));
  }

  ev.put(std::move(seeds));
  ev.put(std::move(sources), "sourceTracks");
}

void HGCTrackSeedEDProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("InputCollection", edm::InputTag("generalTracks"));
  desc.add<std::string>("cut",
                        "1.48 < abs(eta) < 3.0 && pt > 1. && quality(\"highPurity\") && "
                        "hitPattern().numberOfLostHits(\"MISSING_OUTER_HITS\") < 5");
  desc.add<bool>("alongMomentum", true);
  descriptions.add("hgcTrackSeedEDProducer", desc);
}

DEFINE_FWK_MODULE(HGCTrackSeedEDProducer);
