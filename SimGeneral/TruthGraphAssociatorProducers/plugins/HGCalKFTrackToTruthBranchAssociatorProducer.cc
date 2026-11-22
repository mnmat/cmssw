// Original author: Mark Matthewman (CERN) <mark.matthewman@cern.ch>

// Minimal copy of the hit-based path of AllRecoToTruthBranchAssociatorsProducer, producing
// truth-branch association maps for the HGCAL Kalman-filter tracks. The shared producer is
// not usable for them: its reco::Track domain is bound to the tracker channel and to the
// tracker hit adapter, which drops every HGCTrackingRecHit. Integrating HGCAL tracks into
// the shared producer needs further discussion; until then this module writes the same
// products, under the same instance labels, so downstream consumers do not change.
//
// Metric: SharedHits on the Calo channel. One reco hit is one HGCAL cell (DetId).

#include <memory>
#include <span>
#include <string>
#include <vector>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/HGCalReco/interface/HGCTrackingRecHit.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "SimDataFormats/Associations/interface/TICLAssociationMap.h"

#include "PhysicsTools/TruthInfo/interface/Branch.h"
#include "PhysicsTools/TruthInfo/interface/BranchHitAssociator.h"
#include "PhysicsTools/TruthInfo/interface/RecoHitAdapters.h"
#include "SimDataFormats/TruthInfo/interface/Graph.h"
#include "SimDataFormats/TruthInfo/interface/LogicalGraphHitIndex.h"

namespace {
  // reco::Track built by the HGCAL KF -> its valid HGCAL rechits, one unit-weight entry per
  // cell. Missing layers are InvalidTrackingRecHits and are skipped. The result is sorted
  // and has one entry per DetId, as the merge-join in BranchHitAssociator requires.
  std::vector<truth::RecoHit> hgcalRecoHits(reco::Track const& track) {
    std::vector<truth::RecoHit> hits;
    hits.reserve(track.recHitsSize());
    for (auto it = track.recHitsBegin(); it != track.recHitsEnd(); ++it) {
      auto const* hit = dynamic_cast<HGCTrackingRecHit const*>(&**it);
      if (hit == nullptr || !hit->isValid()) {
        continue;
      }
      hits.push_back(truth::RecoHit{.detId = hit->objRef()->id().rawId()});
    }
    truth::sortAndCoalesce(hits);
    return hits;
  }
}  // namespace

class HGCalKFTrackToTruthBranchAssociatorProducer : public edm::global::EDProducer<> {
public:
  explicit HGCalKFTrackToTruthBranchAssociatorProducer(edm::ParameterSet const&);
  void produce(edm::StreamID, edm::Event&, edm::EventSetup const&) const override;
  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  using MapType = ticl::TICLAssociationMap<ticl::mapWithSharedEnergyAndScore>;
  static constexpr auto kMetric = truth::BranchHitAssociator::Metric::SharedHits;
  static constexpr auto kChannel = truth::HitChannel::Calo;

  struct WorkingPoint {
    std::string name;
    float reverseWeight;
    float maxReverseScore;
    bool adaptive;
  };

  const edm::EDGetTokenT<truth::Graph> graphToken_;
  const edm::EDGetTokenT<truth::LogicalGraphHitIndex> hitIndexToken_;
  // The selector-passing candidate roots and the subset a reco object may be assigned to,
  // both from the shared TruthBranchTargetsProducer.
  const edm::EDGetTokenT<std::vector<unsigned int>> targetsToken_;
  const edm::EDGetTokenT<std::vector<unsigned int>> assignableTargetsToken_;
  std::vector<std::pair<std::string, edm::EDGetTokenT<std::vector<reco::Track>>>> recoTokens_;
  std::vector<WorkingPoint> workingPoints_;
};

HGCalKFTrackToTruthBranchAssociatorProducer::HGCalKFTrackToTruthBranchAssociatorProducer(
    edm::ParameterSet const& cfg)
    : graphToken_(consumes<truth::Graph>(cfg.getParameter<edm::InputTag>("src"))),
      hitIndexToken_(consumes<truth::LogicalGraphHitIndex>(cfg.getParameter<edm::InputTag>("hitIndex"))),
      targetsToken_(consumes<std::vector<unsigned int>>(cfg.getParameter<edm::InputTag>("targetsSrc"))),
      assignableTargetsToken_(
          consumes<std::vector<unsigned int>>(cfg.getParameter<edm::InputTag>("assignableTargetsSrc"))) {
  const auto names = cfg.getParameter<std::vector<std::string>>("workingPointNames");
  const auto weights = cfg.getParameter<std::vector<float>>("adaptiveReverseWeight");
  const auto ceilings = cfg.getParameter<std::vector<float>>("adaptiveMaxReverseScore");
  if (names.size() != weights.size() || names.size() != ceilings.size()) {
    throw cms::Exception("Configuration")
        << "workingPointNames, adaptiveReverseWeight and adaptiveMaxReverseScore must have the same length";
  }
  if (names.empty()) {
    throw cms::Exception("Configuration") << "workingPointNames is empty";
  }
  for (std::size_t i = 0; i < names.size(); ++i) {
    // "Fixed" means the plain per-root match; every other point drives the adaptive climb.
    workingPoints_.push_back({names[i], weights[i], ceilings[i], names[i] != "Fixed"});
  }

  for (auto const& tag : cfg.getParameter<std::vector<edm::InputTag>>("recoCollections")) {
    // Same key rule as the shared producer: label and instance joined by an underscore.
    std::string key = tag.label();
    if (!tag.instance().empty()) {
      key += "_" + tag.instance();
    }
    recoTokens_.emplace_back(key, consumes<std::vector<reco::Track>>(tag));
    for (auto const& wp : workingPoints_) {
      produces<MapType>(key + "RecoToTruth" + wp.name);
    }
    produces<MapType>(key + "TruthToReco");
  }
}

void HGCalKFTrackToTruthBranchAssociatorProducer::produce(edm::StreamID,
                                                          edm::Event& event,
                                                          edm::EventSetup const&) const {
  auto const& graph = event.get(graphToken_);
  auto const& hitIndex = event.get(hitIndexToken_);
  auto const& selectedRoots = event.get(targetsToken_);

  const unsigned int nBranches = graph.nParticles();
  std::vector<uint8_t> isAssignable(nBranches, 0);
  for (const unsigned int id : event.get(assignableTargetsToken_)) {
    if (id < nBranches) {
      isAssignable[id] = 1;
    }
  }
  auto assignable = [&](truth::BranchMatch const& m) {
    return m.rootParticleId < nBranches && isAssignable[m.rootParticleId] != 0;
  };

  // Orders a particle before an ancestor that owns the same cells. Declared before the
  // associator, which keeps a span of it.
  const std::vector<uint32_t> generations = truth::particleGenerations(graph);
  const truth::BranchHitAssociator associator(hitIndex,
                                              selectedRoots,
                                              kMetric,
                                              kChannel,
                                              /*emptyRootsMeansAll=*/false,
                                              truth::BranchHitAssociator::kAllDetectors,
                                              /*recHitEnergies=*/nullptr,
                                              generations);

  std::vector<truth::BranchMatch> assignableMatches;
  std::vector<truth::BranchMatch> fixedRow;

  for (auto const& [key, token] : recoTokens_) {
    auto const& tracks = event.get(token);
    const unsigned int nReco = tracks.size();

    std::vector<std::unique_ptr<MapType>> recoToTruthPerWp;
    recoToTruthPerWp.reserve(workingPoints_.size());
    for (std::size_t w = 0; w < workingPoints_.size(); ++w) {
      recoToTruthPerWp.push_back(std::make_unique<MapType>(nReco));
    }
    auto truthToReco = std::make_unique<MapType>(nBranches);

    for (unsigned int i = 0; i < nReco; ++i) {
      const auto hits = hgcalRecoHits(tracks[i]);
      if (hits.empty()) {
        continue;
      }
      const auto matches = associator.bestBranches(std::span<const truth::RecoHit>(hits));

      // An adaptive point may only answer with an assignable root. The fixed point keeps
      // the barred roots too, after every assignable one, so its row [0] is a detector
      // particle whenever one matches. Both keep the ascending-score order.
      assignableMatches.clear();
      for (auto const& match : matches) {
        if (assignable(match)) {
          assignableMatches.push_back(match);
        }
      }
      fixedRow = assignableMatches;
      for (auto const& match : matches) {
        if (!assignable(match)) {
          fixedRow.push_back(match);
        }
      }

      // RECO to TRUTH: score = 1 - reco purity (fraction of the track's hits not covered).
      for (std::size_t w = 0; w < workingPoints_.size(); ++w) {
        auto const& wp = workingPoints_[w];
        if (wp.adaptive) {
          const auto match = truth::BranchHitAssociator::bestAdaptiveBranch(
              std::span<const truth::BranchMatch>(assignableMatches), wp.reverseWeight, wp.maxReverseScore);
          if (match.rootParticleId != truth::BranchMatch::kInvalidRoot) {
            recoToTruthPerWp[w]->insert(i, match.rootParticleId, match.sharedEnergy, match.score);
          }
        } else {
          for (auto const& match : fixedRow) {
            recoToTruthPerWp[w]->insert(i, match.rootParticleId, match.sharedEnergy, match.score);
          }
        }
      }

      // TRUTH to RECO: no adaptive climb; value = shared hit count, score = 1 - truth
      // purity (fraction of the branch's cells not on the track).
      for (auto const& match : matches) {
        if (match.rootParticleId >= nBranches) {
          throw cms::Exception("Configuration")
              << "the hit index and the truth graph disagree on the particle count for '" << key
              << "'. Both products must come from the same event.";
        }
        truthToReco->insert(match.rootParticleId, i, match.sharedEnergy, match.reverseScore);
      }
    }

    // RecoToTruth rows keep their fill order (assignable first, ascending score); sorting
    // them would undo that.
    for (std::size_t w = 0; w < workingPoints_.size(); ++w) {
      event.put(std::move(recoToTruthPerWp[w]), key + "RecoToTruth" + workingPoints_[w].name);
    }
    truthToReco->sort(truth::byAscendingScore);
    event.put(std::move(truthToReco), key + "TruthToReco");
  }
}

void HGCalKFTrackToTruthBranchAssociatorProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("truthLogicalGraphProducer"));
  desc.add<edm::InputTag>("hitIndex", edm::InputTag("truthLogicalGraphHitIndexProducer"));
  desc.add<std::vector<edm::InputTag>>("recoCollections",
                                       {edm::InputTag("ticlTrackstersKalmanFilter", "HGCALTracks")});
  desc.add<edm::InputTag>("targetsSrc", edm::InputTag("truthBranchTargets", "selectedRoots"))
      ->setComment("Candidate roots from TruthBranchTargetsProducer");
  desc.add<edm::InputTag>("assignableTargetsSrc", edm::InputTag("truthBranchTargets", "assignableRoots"))
      ->setComment("The roots an adaptive working point may answer with; the Fixed point also keeps the barred ones");
  desc.add<std::vector<std::string>>("workingPointNames", {"Fixed"})
      ->setComment(
          "One RecoToTruth product per name. \"Fixed\" is the plain per-root match; any other name is adaptive, "
          "driven by its adaptiveReverseWeight and adaptiveMaxReverseScore entries");
  desc.add<std::vector<float>>("adaptiveReverseWeight", {0.f});
  desc.add<std::vector<float>>("adaptiveMaxReverseScore", {0.f});
  descriptions.add("hgcalKFTrackToTruthBranchAssociators", desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(HGCalKFTrackToTruthBranchAssociatorProducer);
