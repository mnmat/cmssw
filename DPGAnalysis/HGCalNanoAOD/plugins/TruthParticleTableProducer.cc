// Original author: Mark Matthewman (CERN) <mark.matthewman@cern.ch>
//
// NanoAOD table of truth-graph particles, for validating reco objects against the
// truth-branch association maps offline.
//
// The per-particle quantities follow Validation/TruthInfo/plugins/TruthBranchRecoValidator.cc
// so that a Python analysis on this table and the DQM validator agree on definitions: the
// calorimeter-entrance eta, the projection onto a level antichain, the branch footprint and
// the production-vertex quantities. What the validator decides from them (Individual,
// Duplicate, Split, Lost, fake) is NOT computed here: that stays in Python so the cuts can
// be changed without reprocessing.
//
// Rows. Writing every particle is not affordable at PU200, so the table holds the union of
//   - the members of every configured truth level (the efficiency denominators),
//   - every particle a configured RecoToTruth association map links to (the candidates a
//     reco object was matched to, so every link index has a row),
//   - optionally every particle with a configured |pdgId| that crossed the calorimeter
//     boundary (e.g. pileup muons, which the signal-only levels do not contain),
//   - optionally the ancestors of all of the above, so the decay chain can be walked.
// Rows are ordered by graph id. Links store graph ids, so join on the graphId column.

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>

#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "PhysicsTools/TruthInfo/interface/SubgraphHitView.h"
#include "SimDataFormats/Associations/interface/TICLAssociationMap.h"
#include "SimDataFormats/TruthInfo/interface/Graph.h"
#include "SimDataFormats/TruthInfo/interface/LogicalGraphHitIndex.h"
#include "SimDataFormats/TruthInfo/interface/Particle.h"
#include "SimDataFormats/TruthInfo/interface/VertexData.h"

namespace {
  // Same sentinel as truth::kNoCaloEntry in Validation/TruthInfo, for a quantity that has
  // no value (no direction, no production vertex, no calorimeter entry).
  constexpr float kNoValue = -999.f;
  constexpr const char* kLevelInstance = "truthToRecoTargets";

  std::string capitalized(std::string name) {
    if (!name.empty()) {
      name[0] = std::toupper(static_cast<unsigned char>(name[0]));
    }
    return name;
  }

  // Membership mask of a sorted-or-not id list, sized to the graph.
  std::vector<uint8_t> maskOf(std::vector<unsigned int> const& ids, uint32_t n) {
    std::vector<uint8_t> mask(n, 0);
    for (const unsigned int id : ids) {
      if (id < n) {
        mask[id] = 1;
      }
    }
    return mask;
  }
}  // namespace

class TruthParticleTableProducer : public edm::global::EDProducer<> {
public:
  explicit TruthParticleTableProducer(edm::ParameterSet const&);
  void produce(edm::StreamID, edm::Event&, edm::EventSetup const&) const override;
  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  using MapType = ticl::TICLAssociationMap<ticl::mapWithSharedEnergyAndScore>;

  struct LevelEntry {
    std::string name;
    edm::EDGetTokenT<std::vector<unsigned int>> targetsToken;
    edm::EDGetTokenT<std::vector<unsigned int>> eligibilityToken;
  };

  const std::string name_;
  const std::string doc_;
  const edm::EDGetTokenT<truth::Graph> graphToken_;
  const edm::EDGetTokenT<truth::LogicalGraphHitIndex> hitIndexToken_;
  const edm::EDGetTokenT<std::vector<unsigned int>> selectedToken_;
  const edm::EDGetTokenT<std::vector<unsigned int>> assignableToken_;
  const edm::EDGetTokenT<std::vector<unsigned int>> denominatorToken_;
  const edm::EDGetTokenT<std::vector<unsigned int>> signalSeedsToken_;
  std::vector<LevelEntry> levels_;
  std::size_t rootLevelIndex_ = 0;
  std::vector<edm::EDGetTokenT<MapType>> recoToTruthTokens_;
  std::vector<int> keepPdgIds_;
  const bool keepAncestors_;
};

TruthParticleTableProducer::TruthParticleTableProducer(edm::ParameterSet const& cfg)
    : name_(cfg.getParameter<std::string>("name")),
      doc_(cfg.getParameter<std::string>("doc")),
      graphToken_(consumes<truth::Graph>(cfg.getParameter<edm::InputTag>("src"))),
      hitIndexToken_(consumes<truth::LogicalGraphHitIndex>(cfg.getParameter<edm::InputTag>("hitIndex"))),
      selectedToken_(consumes<std::vector<unsigned int>>(
          edm::InputTag(cfg.getParameter<std::string>("targetsProducer"), "selectedRoots"))),
      assignableToken_(consumes<std::vector<unsigned int>>(
          edm::InputTag(cfg.getParameter<std::string>("targetsProducer"), "assignableRoots"))),
      denominatorToken_(consumes<std::vector<unsigned int>>(
          edm::InputTag(cfg.getParameter<std::string>("targetsProducer"), "denominatorRoots"))),
      signalSeedsToken_(consumes<std::vector<unsigned int>>(
          edm::InputTag(cfg.getParameter<std::string>("targetsProducer"), "signalSeeds"))),
      keepAncestors_(cfg.getParameter<bool>("keepAncestors")) {
  const auto targetsProducer = cfg.getParameter<std::string>("targetsProducer");
  const auto levelNames = cfg.getParameter<std::vector<std::string>>("truthLevels");
  const auto rootLevel = cfg.getParameter<std::string>("rootLevel");
  if (levelNames.empty()) {
    throw cms::Exception("Configuration") << "TruthParticleTableProducer: truthLevels is empty";
  }
  // Must match the levels TruthBranchTargetsProducer was configured with, or the products
  // do not exist and the job fails on the first event.
  for (auto const& level : levelNames) {
    const std::string instance = std::string(kLevelInstance) + capitalized(level);
    levels_.push_back({level,
                       consumes<std::vector<unsigned int>>(edm::InputTag(targetsProducer, instance)),
                       consumes<std::vector<unsigned int>>(edm::InputTag(targetsProducer, instance + "Eligibility"))});
  }
  auto const it = std::find(levelNames.begin(), levelNames.end(), rootLevel);
  if (it == levelNames.end()) {
    throw cms::Exception("Configuration")
        << "TruthParticleTableProducer: rootLevel '" << rootLevel << "' is not in truthLevels";
  }
  rootLevelIndex_ = std::distance(levelNames.begin(), it);

  for (auto const& tag : cfg.getParameter<std::vector<edm::InputTag>>("recoToTruth")) {
    recoToTruthTokens_.push_back(consumes<MapType>(tag));
  }
  for (const int pdgId : cfg.getParameter<std::vector<int>>("keepPdgIds")) {
    keepPdgIds_.push_back(std::abs(pdgId));
  }
  produces<nanoaod::FlatTable>();
}

void TruthParticleTableProducer::produce(edm::StreamID, edm::Event& event, edm::EventSetup const&) const {
  auto const& graph = event.get(graphToken_);
  auto const& hitIndexProduct = event.get(hitIndexToken_);
  // Per event: the view caches coalesced subgraphs and is not thread safe.
  truth::SubgraphHitView hitIndex(hitIndexProduct);
  const uint32_t nParticles = graph.nParticles();

  const auto isSelected = maskOf(event.get(selectedToken_), nParticles);
  const auto isAssignable = maskOf(event.get(assignableToken_), nParticles);
  const auto isDenominator = maskOf(event.get(denominatorToken_), nParticles);
  const auto isSignalSeed = maskOf(event.get(signalSeedsToken_), nParticles);

  // Level membership and eligibility, one pair of columns per level.
  std::vector<std::vector<uint8_t>> isLevel(levels_.size(), std::vector<uint8_t>(nParticles, 0));
  std::vector<std::vector<uint32_t>> eligibility(levels_.size(), std::vector<uint32_t>(nParticles, 0));
  std::vector<uint8_t> keep(nParticles, 0);
  for (std::size_t l = 0; l < levels_.size(); ++l) {
    auto const& targets = event.get(levels_[l].targetsToken);
    auto const& elig = event.get(levels_[l].eligibilityToken);
    for (std::size_t t = 0; t < targets.size(); ++t) {
      const unsigned int id = targets[t];
      if (id >= nParticles) {
        continue;
      }
      isLevel[l][id] = 1;
      eligibility[l][id] = (t < elig.size()) ? elig[t] : 0u;
      keep[id] = 1;
    }
  }

  // Every particle a reco object links to, so that no link index is left without a row.
  for (auto const& token : recoToTruthTokens_) {
    for (auto const& row : event.get(token).getMap()) {
      for (auto const& element : row) {
        if (element.index() < nParticles) {
          keep[element.index()] = 1;
        }
      }
    }
  }

  // Calorimeter entrance, as in TruthBranchRecoValidator: the boundary crossing
  // (checkpoint 0) of the most energetic particle of the branch to reach the calorimeter,
  // propagated up the production edges so a root inherits it from its descendants.
  // Back-scattered tracks crossed inward and are skipped.
  std::vector<float> caloEntryEta(nParticles, kNoValue);
  std::vector<uint8_t> crossesBoundary(nParticles, 0);
  {
    std::vector<double> caloEntryEnergy(nParticles, -1.);
    std::vector<uint32_t> worklist;
    worklist.reserve(nParticles);
    for (uint32_t p = 0; p < nParticles; ++p) {
      if (graph.particles()[p].backscattered) {
        continue;
      }
      const auto crossing = truth::Particle(&graph, p).checkpoint(0);
      if (!crossing.has_value()) {
        continue;
      }
      crossesBoundary[p] = 1;
      caloEntryEnergy[p] = crossing->momentum.energy();
      caloEntryEta[p] = crossing->position.eta();
      worklist.push_back(p);
    }
    for (std::size_t head = 0; head < worklist.size(); ++head) {
      const uint32_t p = worklist[head];
      for (const uint32_t vertexId : graph.productionVertices(p)) {
        if (vertexId >= graph.nVertices()) {
          continue;
        }
        for (const uint32_t parent : graph.incomingParticles(vertexId)) {
          if (parent < nParticles && caloEntryEnergy[p] > caloEntryEnergy[parent]) {
            caloEntryEnergy[parent] = caloEntryEnergy[p];
            caloEntryEta[parent] = caloEntryEta[p];
            worklist.push_back(parent);
          }
        }
      }
    }
  }

  // Extra species that crossed into the calorimeter, whichever interaction they came from.
  if (!keepPdgIds_.empty()) {
    for (uint32_t p = 0; p < nParticles; ++p) {
      if (crossesBoundary[p] &&
          std::find(keepPdgIds_.begin(), keepPdgIds_.end(), std::abs(graph.particles()[p].pdgId)) !=
              keepPdgIds_.end()) {
        keep[p] = 1;
      }
    }
  }

  // Projection onto the root level, as the validator's dominance projection: for every
  // particle the level member it descends from, or itself if it is one. First writer wins,
  // which also bounds the walk on a graph with a cycle.
  std::vector<int32_t> levelRoot(nParticles, -1);
  {
    std::vector<uint32_t> worklist;
    worklist.reserve(nParticles);
    for (uint32_t id = 0; id < nParticles; ++id) {
      if (isLevel[rootLevelIndex_][id]) {
        levelRoot[id] = static_cast<int32_t>(id);
        worklist.push_back(id);
      }
    }
    for (std::size_t head = 0; head < worklist.size(); ++head) {
      const uint32_t p = worklist[head];
      for (const uint32_t vertexId : graph.decayVertices(p)) {
        if (vertexId >= graph.nVertices()) {
          continue;
        }
        for (const uint32_t child : graph.outgoingParticles(vertexId)) {
          if (child < nParticles && levelRoot[child] < 0) {
            levelRoot[child] = levelRoot[p];
            worklist.push_back(child);
          }
        }
      }
    }
  }

  // Close the row set under ancestry, so parentGraphId always points at a row.
  if (keepAncestors_) {
    std::vector<uint32_t> worklist;
    for (uint32_t id = 0; id < nParticles; ++id) {
      if (keep[id]) {
        worklist.push_back(id);
      }
    }
    for (std::size_t head = 0; head < worklist.size(); ++head) {
      const uint32_t p = worklist[head];
      for (const uint32_t vertexId : graph.productionVertices(p)) {
        if (vertexId >= graph.nVertices()) {
          continue;
        }
        for (const uint32_t parent : graph.incomingParticles(vertexId)) {
          if (parent < nParticles && !keep[parent]) {
            keep[parent] = 1;
            worklist.push_back(parent);
          }
        }
      }
    }
  }

  std::vector<uint32_t> rows;
  for (uint32_t id = 0; id < nParticles; ++id) {
    if (keep[id]) {
      rows.push_back(id);
    }
  }
  const std::size_t n = rows.size();

  std::vector<uint32_t> graphId(n), nCaloCells(n), nDirectCaloCells(n), nMuonHits(n);
  std::vector<int32_t> pdgId(n), parentGraphId(n), rootGraphId(n);
  std::vector<int16_t> status(n), depth(n);
  std::vector<float> charge(n), pt(n), eta(n), phi(n), energy(n), caloEta(n);
  std::vector<float> vx(n), vy(n), vz(n), vertpos(n), dxy(n), dz(n);
  std::vector<uint8_t> vertexReason(n);
  std::vector<bool> isSignal(n), backscattered(n), synthetic(n), crosses(n);
  std::vector<bool> selected(n), assignable(n), denominator(n), signalSeed(n);
  std::vector<std::vector<bool>> levelCol(levels_.size(), std::vector<bool>(n));
  std::vector<std::vector<uint32_t>> eligCol(levels_.size(), std::vector<uint32_t>(n));

  for (std::size_t r = 0; r < n; ++r) {
    const uint32_t id = rows[r];
    auto const& data = graph.particles()[id];
    const truth::Particle particle(&graph, id);
    auto const& p4 = data.momentum;
    // A resonance in its pre-ISR copy has exactly zero pt and so no direction.
    const bool hasDirection = p4.pt() > 0.;

    graphId[r] = id;
    pdgId[r] = data.pdgId;
    status[r] = data.status;
    charge[r] = particle.charge();
    pt[r] = p4.pt();
    eta[r] = hasDirection ? p4.eta() : kNoValue;
    phi[r] = hasDirection ? p4.phi() : kNoValue;
    energy[r] = p4.energy();
    caloEta[r] = caloEntryEta[id];
    isSignal[r] = data.isSignal();
    backscattered[r] = data.backscattered;
    synthetic[r] = data.isSynthetic();
    crosses[r] = crossesBoundary[id] != 0;
    depth[r] = static_cast<int16_t>(std::min<uint32_t>(particle.ancestorCount(), 32767u));

    // The branch footprint: the truth analogue of a reco object's hit count. nCaloCells
    // is the denominator of the shared-hits truth purity (coverage = nShared / nCaloCells).
    nCaloCells[r] = hitIndex.subgraphHits(truth::HitChannel::Calo, id).size();
    nDirectCaloCells[r] = hitIndex.directHits(truth::HitChannel::Calo, id).size();
    nMuonHits[r] = hitIndex.subgraphHits(truth::HitChannel::Muon, id).size();

    // First parent through the first production vertex; -1 for a beam particle.
    parentGraphId[r] = -1;
    vx[r] = vy[r] = vz[r] = vertpos[r] = dxy[r] = dz[r] = kNoValue;
    vertexReason[r] = 0;
    auto const production = graph.productionVertices(id);
    if (!production.empty() && production.front() < graph.nVertices()) {
      const uint32_t vertexId = production.front();
      auto const parents = graph.incomingParticles(vertexId);
      if (!parents.empty()) {
        parentGraphId[r] = static_cast<int32_t>(parents.front());
      }
      auto const& vdata = graph.vertices()[vertexId];
      auto const& pos = vdata.position;
      vx[r] = pos.x();
      vy[r] = pos.y();
      vz[r] = pos.z();
      vertpos[r] = std::sqrt(pos.x() * pos.x() + pos.y() * pos.y());
      vertexReason[r] = vdata.reason;
      // Transverse and longitudinal impact parameter of the branch direction with respect
      // to the origin, the truth counterpart of the track dxy and dz.
      if (hasDirection) {
        dxy[r] = (-pos.x() * p4.py() + pos.y() * p4.px()) / p4.pt();
        dz[r] = pos.z() - (pos.x() * p4.px() + pos.y() * p4.py()) / p4.pt() * (p4.pz() / p4.pt());
      }
    }

    rootGraphId[r] = levelRoot[id];
    selected[r] = isSelected[id] != 0;
    assignable[r] = isAssignable[id] != 0;
    denominator[r] = isDenominator[id] != 0;
    signalSeed[r] = isSignalSeed[id] != 0;
    for (std::size_t l = 0; l < levels_.size(); ++l) {
      levelCol[l][r] = isLevel[l][id] != 0;
      eligCol[l][r] = eligibility[l][id];
    }
  }

  auto table = std::make_unique<nanoaod::FlatTable>(n, name_, false);
  table->setDoc(doc_);
  table->addColumn<uint32_t>("graphId", graphId, "truth graph particle id; association links index this");
  table->addColumn<int32_t>("pdgId", pdgId, "PDG id");
  table->addColumn<int16_t>("status", status, "generator status (0 for SIM-only particles)");
  table->addColumn<float>("charge", charge, "electric charge");
  table->addColumn<float>("pt", pt, "pt of the nominal momentum (GEN if available, else SimTrack)", 12);
  table->addColumn<float>("eta", eta, "eta; -999 when pt is 0", 12);
  table->addColumn<float>("phi", phi, "phi; -999 when pt is 0", 12);
  table->addColumn<float>("energy", energy, "energy", 12);
  table->addColumn<float>(
      "caloEta", caloEta, "eta at the calorimeter entrance of the most energetic crossing in the branch; -999 if none", 12);
  table->addColumn<bool>("crossesCaloBoundary", crosses, "this particle itself crossed the tracker-calo boundary outward");
  table->addColumn<bool>("isSignal", isSignal, "from the signal interaction (false = pileup)");
  table->addColumn<bool>("backscattered", backscattered, "Geant4 back-scattered (crossed the boundary inward)");
  table->addColumn<bool>("isSynthetic", synthetic, "a node the graph invented (connector, signal stand-in)");
  table->addColumn<int16_t>("depth", depth, "number of ancestors");
  table->addColumn<int32_t>("parentGraphId", parentGraphId, "graph id of the first parent; -1 if none");
  table->addColumn<float>("vx", vx, "production vertex x [cm]; -999 if none", 12);
  table->addColumn<float>("vy", vy, "production vertex y [cm]; -999 if none", 12);
  table->addColumn<float>("vz", vz, "production vertex z [cm]; -999 if none", 12);
  table->addColumn<float>("vertpos", vertpos, "production vertex transverse radius [cm]; -999 if none", 12);
  table->addColumn<float>("dxy", dxy, "transverse impact parameter w.r.t. the origin [cm]; -999 if undefined", 12);
  table->addColumn<float>("dz", dz, "longitudinal impact parameter w.r.t. the origin [cm]; -999 if undefined", 12);
  table->addColumn<uint8_t>("vertexReason", vertexReason, "truth::VertexReason of the production vertex");
  table->addColumn<uint32_t>(
      "nCaloCells", nCaloCells, "Calo-channel cells (unique DetIds) of the branch: the particle and its descendants");
  table->addColumn<uint32_t>("nDirectCaloCells", nDirectCaloCells, "Calo-channel cells of the particle itself");
  table->addColumn<uint32_t>("nMuonHits", nMuonHits, "Muon-channel entries of the branch (chamber/roll level)");
  table->addColumn<bool>("isSelected", selected, "in truthBranchTargets selectedRoots: a reco object can match it");
  table->addColumn<bool>("isAssignable", assignable, "in assignableRoots: a reco object may be assigned to it");
  table->addColumn<bool>("isDenominator", denominator, "in denominatorRoots: passes every selector cut");
  table->addColumn<bool>("isSignalSeed", signalSeed, "in signalSeeds");
  for (std::size_t l = 0; l < levels_.size(); ++l) {
    auto const& level = levels_[l].name;
    table->addColumn<bool>("isLevel_" + level, levelCol[l], "member of the " + level + " truth level (antichain)");
    table->addColumn<uint32_t>("eligibility_" + level,
                               eligCol[l],
                               "for " + level + " members: the kinematic selector cut it fails (0 = passes all); "
                               "0 for non-members");
  }
  table->addColumn<int32_t>("rootGraphId",
                            rootGraphId,
                            "graph id of the " + levels_[rootLevelIndex_].name +
                                " member this particle is or descends from; -1 if none (above the level, or "
                                "in no member's branch)");
  event.put(std::move(table));
}

void TruthParticleTableProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::string>("name", "TruthParticle");
  desc.add<std::string>("doc", "truth graph particles: level members, association candidates and their ancestors");
  desc.add<edm::InputTag>("src", edm::InputTag("truthLogicalGraphProducer"));
  desc.add<edm::InputTag>("hitIndex", edm::InputTag("truthLogicalGraphHitIndexProducer"));
  desc.add<std::string>("targetsProducer", "truthBranchTargets");
  desc.add<std::vector<std::string>>("truthLevels", {"caloBoundary"})
      ->setComment("Must be the truthLevels the targets producer was configured with (or a subset)");
  desc.add<std::string>("rootLevel", "caloBoundary")
      ->setComment("The level rootGraphId projects onto; the validator's dominanceLevel");
  desc.add<std::vector<edm::InputTag>>("recoToTruth", {})
      ->setComment("RecoToTruth association maps whose linked particles get a row; use the Fixed working point");
  desc.add<std::vector<int>>("keepPdgIds", {})
      ->setComment("Also keep every particle with this |pdgId| that crossed the calorimeter boundary, pileup included");
  desc.add<bool>("keepAncestors", true)->setComment("Also keep every ancestor of a kept particle");
  descriptions.add("truthParticleTable", desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TruthParticleTableProducer);
