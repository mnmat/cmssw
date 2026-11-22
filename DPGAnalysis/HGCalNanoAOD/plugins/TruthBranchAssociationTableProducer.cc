// Original author: Mark Matthewman (CERN) <mark.matthewman@cern.ch>
//
// NanoAOD link tables for the truth-branch association maps
// (SimGeneral/TruthGraphAssociatorProducers). These maps carry no collection types, so the
// existing TracksterTracksterEnergyScoreFlatTableProducer cannot read them. In a
// RecoToTruth map the main table has one row per reco object and each link's index is a
// truth-graph particle id (join on TruthParticle_graphId).

#include "PhysicsTools/NanoAOD/interface/AssociationMapFlatTableProducer.h"
#include "SimDataFormats/Associations/interface/TICLAssociationMap.h"

typedef AssociationOneToManyFlatTableProducer<ticl::TICLAssociationMap<ticl::mapWithSharedEnergyAndScore>>
    TruthBranchAssociationFlatTableProducer;

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TruthBranchAssociationFlatTableProducer);
