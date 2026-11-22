import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *

from SimTracker.TrackAssociation.trackingParticleRecoTrackAssociationTables_cff import (
    trackingParticleRecoTrackAssociationTable,
    recoTrackTrackingParticleAssociationTable,
)

tpTable = cms.EDProducer("SimpleTrackingParticleFlatTableProducer",
    src = cms.InputTag("mix", "MergedTrackTruth"),
    cut = cms.string(""),              # keep empty: row == TP index used by TPAssoc tables
    name = cms.string("TP"), doc = cms.string("TrackingParticles"),
    singleton = cms.bool(False), extension = cms.bool(False),
    skipNonExistingSrc = cms.bool(False),
    variables = cms.PSet(
        pt = Var("pt", "float"), eta = Var("eta", "float"), phi = Var("phi", "float"),
        pdgId = Var("pdgId", "int"), charge = Var("charge", "int"),
        g4TrackId = Var("g4Tracks()[0].trackId()", "uint", doc="first SimTrack id"),
        eventId   = Var("eventId().rawId()", "uint", doc="EncodedEventId raw"),
    ),
)

cpTable = cms.EDProducer("SimpleCaloParticleFlatTableProducer",
    src = cms.InputTag("mix", "MergedCaloTruth"),
    cut = cms.string(""),              # keep empty: row == CP index used by SimTrackster/assoc tables
    name = cms.string("CP"), doc = cms.string("CaloParticles"),
    singleton = cms.bool(False), extension = cms.bool(False),
    skipNonExistingSrc = cms.bool(False),
    variables = cms.PSet(
        energy = Var("energy", "float"), pt = Var("pt", "float"),
        eta = Var("eta", "float"), phi = Var("phi", "float"),
        pdgId = Var("pdgId", "int"), simEnergy = Var("simEnergy", "float"),
        g4TrackId = Var("g4Tracks()[0].trackId()", "uint"),
        eventId   = Var("g4Tracks()[0].eventId().rawId()", "uint"),
    ),
)


# HGC Tracks table for NanoAOD 
hgcTracksTable = cms.EDProducer(
    "SimpleTrackFlatTableProducer",
    skipNonExistingSrc=cms.bool(True),
    src=cms.InputTag("ticlTrackstersKalmanFilter","HGCALTracks","RECO"),
    cut=cms.string(""),
    name=cms.string("HGCTrack"),
    doc=cms.string("HGCTracks"),
    extension=cms.bool(False),
    variables=cms.PSet(
        # Basic kinematics
        pt=Var("pt()", "float", doc="track p_T [GeV]"),
        p=Var("p()", "float", doc="track momentum magnitude [GeV]"),
        eta=Var("eta()", "float", doc="track pseudorapidity"),
        phi=Var("phi()", "float", doc="track phi angle [rad]"),
        charge=Var("charge()", "int", doc="track charge"),
        trackLambda=Var("lambda()", "float", doc="track lambda"),
        
        # Position 
        vx=Var("vx()", "float", doc="track vertex x [cm]"),
        vy=Var("vy()", "float", doc="track vertex y [cm]"),
        vz=Var("vz()", "float", doc="track vertex z [cm]"),
        
        # Quality metrics
        nhits=Var("numberOfValidHits()", "uint16", doc="number of valid hits"),
        missingOuterHits=Var("missingOuterHits()", "uint8", doc="number of missing outer hits"),
        
        # Error parameters 
        ptErr=Var("ptError()", "float", doc="track p_T error [GeV]"),
        etaErr=Var("etaError()", "float", doc="track eta error"),
        phiErr=Var("phiError()", "float", doc="track phi error"),
        lambdaErr=Var("lambdaError()", "float", doc="track lambda error"),
        qoverpErr=Var("qoverpError()", "float", doc="q/p error"),
    ),
    externalVariables = cms.PSet(
        seedingTrackIdx = ExtVar(cms.InputTag("ticlTrackstersKalmanFilter","HGCSeedingTrackIds","RECO"),"int", doc="index of the seeding track in the input collection"),
    )
    # externalVariables = cms.PSet(
    #     isMuon = ExtVar(cms.InputTag("trackMuonInfo", "isMuon"),"int", doc="1 if PFMuonAlgo::isMuon, 0 otherwise, -1 if no muon"),
    #     isTrackerMuon = ExtVar(cms.InputTag("trackMuonInfo", "isTrackerMuon"),"int", doc="1 if tracker muon, -1 if no muon"),
    #     muon_dt_hits = ExtVar(cms.InputTag("trackMuonInfo", "muonDtHits"),"int", doc="ST muon DT hits (-1 if no muon)"),
    #     muon_csc_hits = ExtVar(cms.InputTag("trackMuonInfo", "muonCscHits"), "int", doc="ST muon CSC hits (-1 if no muon)"),
    #     muon_type = ExtVar(cms.InputTag("trackMuonInfo", "muonType"), "int", doc="muon type (-1 if no muon)"),
    # )
)

hgcTrackRecHitTableProducer = cms.EDProducer(
    "HGCTrackRecHitTableProducer",
    skipNonExistingSrc=cms.bool(True),
    hgcTracks=cms.InputTag("ticlTrackstersKalmanFilter","HGCALTracks","RECO"),
    hitTableName=cms.string("HGCTrackRecHit"),
    trackTableName=cms.string("HGCTrackRecHitLinks"),
    precision=cms.int32(6))

cpRecHitTableProducer = cms.EDProducer(
    "CPRecHitTableProducer",
    skipNonExistingSrc=cms.bool(True),
    caloParticles=cms.InputTag("mix", "MergedCaloTruth"),
    simVertices = cms.InputTag("g4SimHits"),
    hitTableName=cms.string("CPRecHit"),
    trackTableName=cms.string("CPRecHitLinks"),
    precision=cms.int32(6))

# Sequence for general tracks
hgcalTracksTableSequence = cms.Sequence(
    hgcTracksTable
    + trackingParticleRecoTrackAssociationTable
    + recoTrackTrackingParticleAssociationTable
    + tpTable
    + cpTable
    + hgcTrackRecHitTableProducer
    + cpRecHitTableProducer
)