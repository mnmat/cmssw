import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import Var

# Truth-branch validation of the HGCAL KF tracks in NanoAOD:
#   HGCKFTrack<->truth association (SharedHits, Calo channel)
#   -> link table per track collection (index = truth graph id)
#   -> TruthParticle table (rows = level members + linked particles + ancestors)
# Efficiency, fake rate, duplicates and pollution are computed offline from these tables.
#
# Needs truthLogicalGraphProducer and truthLogicalGraphHitIndexProducer in the job.

from SimGeneral.TruthGraphAssociatorProducers.truthGraphAssociators_cff import truthBranchTargets as _truthBranchTargets

_kfCollections = [
    cms.InputTag("ticlTrackstersKalmanFilter", "HGCALTracks"),
    cms.InputTag("ticlTrackstersKalmanFilterGE0", "HGCALTracks"),
]
_truthLevels = ["caloBoundary", "reconstructableFromSignal"]

def _key(tag):
    return tag.getModuleLabel() + ("_" + tag.getProductInstanceLabel() if tag.getProductInstanceLabel() else "")

hgcKFTruthBranchTargets = _truthBranchTargets.clone(
    truthLevels = _truthLevels,
)

from SimGeneral.TruthGraphAssociatorProducers.hgcalKFTrackToTruthBranchAssociators_cfi import hgcalKFTrackToTruthBranchAssociators as _assoc
hgcalKFTrackToTruthBranchAssociators = _assoc.clone(
    recoCollections = _kfCollections,
    targetsSrc = ("hgcKFTruthBranchTargets", "selectedRoots"),
    assignableTargetsSrc = ("hgcKFTruthBranchTargets", "assignableRoots"),
    workingPointNames = ["Fixed"],
    adaptiveReverseWeight = [0.],
    adaptiveMaxReverseScore = [0.],
)

def _linkTable(tag, name):
    key = _key(tag)
    return cms.EDProducer(
        "TruthBranchAssociationFlatTableProducer",
        src = cms.InputTag("hgcalKFTrackToTruthBranchAssociators", key + "RecoToTruthFixed"),
        name = cms.string(name),
        doc = cms.string(f"{key} -> truth graph particles (Fixed WP): assignable roots first, ascending score"),
        extension = cms.bool(False),
        skipNonExistingSrc = cms.bool(False),
        variables = cms.PSet(),
        collectionVariables = cms.PSet(
            links = cms.PSet(
                name = cms.string(name + "Links"),
                doc = cms.string("association links"),
                useCount = cms.bool(True),
                useOffset = cms.bool(True),
                variables = cms.PSet(
                    index = Var("index", "uint", doc="truth graph id (join on TruthParticle_graphId)"),
                    sharedEnergy = Var("sharedEnergy", "float", doc="number of shared HGCAL cells"),
                    score = Var("score", "float", doc="1 - track purity (fraction of the track's cells not from this branch)"),
                ),
            ),
        ),
    )

hgcKFTrack2TruthTable = _linkTable(_kfCollections[0], "HGCTrack2Truth")
hgcKFTrackGE02TruthTable = _linkTable(_kfCollections[1], "HGCTrackGE02Truth")

from DPGAnalysis.HGCalNanoAOD.truthParticleTable_cfi import truthParticleTable as _truthParticleTable
hgcKFTruthParticleTable = _truthParticleTable.clone(
    targetsProducer = "hgcKFTruthBranchTargets",
    truthLevels = _truthLevels,
    rootLevel = "caloBoundary",
    recoToTruth = [cms.InputTag("hgcalKFTrackToTruthBranchAssociators", _key(t) + "RecoToTruthFixed") for t in _kfCollections],
    keepPdgIds = [13],        # pileup muons too, which the signal-only levels do not hold
    keepAncestors = True,
)

hgcalKFTruthTablesSequence = cms.Sequence(
    hgcKFTruthBranchTargets
    + hgcalKFTrackToTruthBranchAssociators
    + hgcKFTrack2TruthTable
    + hgcKFTrackGE02TruthTable
    + hgcKFTruthParticleTable
)
