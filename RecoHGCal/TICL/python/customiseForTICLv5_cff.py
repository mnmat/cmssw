import FWCore.ParameterSet.Config as cms

from RecoHGCal.TICL.ticlLayerTileProducer_cfi import ticlLayerTileProducer

from RecoHGCal.TICL.CLUE3DEM_cff import *
from RecoHGCal.TICL.CLUE3DHAD_cff import *
from RecoHGCal.TICL.pfTICLProducerV5_cfi import pfTICLProducerV5 as _pfTICLProducerV5

from RecoHGCal.TICL.ticlLayerTileProducer_cfi import ticlLayerTileProducer
from RecoHGCal.TICL.tracksterSelectionTf_cfi import *

from RecoHGCal.TICL.tracksterLinksProducer_cfi import tracksterLinksProducer as _tracksterLinksProducer
from RecoHGCal.TICL.ticlCandidateProducer_cfi import ticlCandidateProducer as _ticlCandidateProducer
from RecoHGCal.Configuration.RecoHGCal_EventContent_cff import customiseForTICLv5EventContent
from RecoHGCal.TICL.iterativeTICL_cff import ticlIterLabels, ticlIterLabelsMerge
from RecoHGCal.TICL.ticlDumper_cfi import ticlDumper
from RecoHGCal.TICL.mergedTrackstersProducer_cfi import mergedTrackstersProducer as _mergedTrackstersProducer
from SimCalorimetry.HGCalAssociatorProducers.TSToSimTSAssociation_cfi import tracksterSimTracksterAssociationLinkingbyCLUE3D as _tracksterSimTracksterAssociationLinkingbyCLUE3D
from SimCalorimetry.HGCalAssociatorProducers.TSToSimTSAssociation_cfi import tracksterSimTracksterAssociationPRbyCLUE3D  as _tracksterSimTracksterAssociationPRbyCLUE3D 
from Validation.HGCalValidation.HGCalValidator_cff import hgcalValidatorv5 
from RecoLocalCalo.HGCalRecProducers.HGCalUncalibRecHit_cfi import HGCalUncalibRecHit
from RecoHGCal.TICL.SimTracksters_cff import ticlSimTracksters

from RecoHGCal.TICL.FastJetStep_cff import ticlTrackstersFastJet
from RecoHGCal.TICL.EMStep_cff import ticlTrackstersEM, ticlTrackstersHFNoseEM
from RecoHGCal.TICL.TrkStep_cff import ticlTrackstersTrk, ticlTrackstersHFNoseTrk
from RecoHGCal.TICL.MIPStep_cff import ticlTrackstersMIP, ticlTrackstersHFNoseMIP
from RecoHGCal.TICL.HADStep_cff import ticlTrackstersHAD, ticlTrackstersHFNoseHAD
from RecoHGCal.TICL.CLUE3DEM_cff import ticlTrackstersCLUE3DEM
from RecoHGCal.TICL.CLUE3DHAD_cff import ticlTrackstersCLUE3DHAD
from RecoHGCal.TICL.CLUE3DHighStep_cff import ticlTrackstersCLUE3DHigh
from RecoHGCal.TICL.TrkEMStep_cff import ticlTrackstersTrkEM, filteredLayerClustersHFNoseTrkEM

from RecoHGCal.TICL.mtdSoAProducer_cfi import mtdSoAProducer as _mtdSoAProducer

def customiseForTICLv5(process, enableDumper = False):

    process.HGCalUncalibRecHit.computeLocalTime = cms.bool(True)
    process.ticlSimTracksters.computeLocalTime = cms.bool(True)

    process.ticlTrackstersFastJet.pluginPatternRecognitionByFastJet.computeLocalTime = cms.bool(True)

    process.ticlTrackstersEM.pluginPatternRecognitionByCA.computeLocalTime = cms.bool(True)
    process.ticlTrackstersHFNoseEM.pluginPatternRecognitionByCA.computeLocalTime = cms.bool(True)

    process.ticlTrackstersTrk.pluginPatternRecognitionByCA.computeLocalTime = cms.bool(True)
    process.ticlTrackstersHFNoseTrk.pluginPatternRecognitionByCA.computeLocalTime = cms.bool(True)

    process.ticlTrackstersMIP.pluginPatternRecognitionByCA.computeLocalTime = cms.bool(True)
    process.ticlTrackstersHFNoseMIP.pluginPatternRecognitionByCA.computeLocalTime = cms.bool(True)

    process.ticlTrackstersHAD.pluginPatternRecognitionByCA.computeLocalTime = cms.bool(True)
    process.ticlTrackstersHFNoseHAD.pluginPatternRecognitionByCA.computeLocalTime = cms.bool(True)

    process.ticlTrackstersCLUE3DHAD.pluginPatternRecognitionByCLUE3D.computeLocalTime = cms.bool(True)
    process.ticlTrackstersCLUE3DEM.pluginPatternRecognitionByCLUE3D.computeLocalTime = cms.bool(True)
    process.ticlTrackstersCLUE3DHigh.pluginPatternRecognitionByCLUE3D.computeLocalTime = cms.bool(True)

    process.ticlTrackstersTrkEM.pluginPatternRecognitionByCA.computeLocalTime = cms.bool(True)
    process.ticlTrackstersHFNoseTrkEM.pluginPatternRecognitionByCA.computeLocalTime = cms.bool(True)

    process.ticlLayerTileTask = cms.Task(ticlLayerTileProducer)

    process.filteredLayerClustersCLUE3DEM.max_layerId = 47
    process.ticlTrackstersCLUE3DEM.pluginPatternRecognitionByCLUE3D = cms.PSet(
        algo_verbosity = cms.int32(0),
        criticalDensity = cms.vdouble(0.9, 0.9, 0.613107466),
        criticalEtaPhiDistance = cms.vdouble(0.025, 0.025, 0.025),
        criticalSelfDensity = cms.vdouble(0.1, 0.355379138, 0.1),
        criticalXYDistance = cms.vdouble(1.0, 2.5, 2.5),
        criticalZDistanceLyr = cms.vint32(7, 7, 4),
        cutHadProb = cms.double(0.5),
        densityEtaPhiDistanceSqr = cms.vdouble(0.0008, 0.0008, 0.0008),
        densityOnSameLayer = cms.bool(False),
        densitySiblingLayers = cms.vint32(5, 2, 5),
        densityXYDistanceSqr = cms.vdouble(4.16576578, 5.0, 5.0),
        doPidCut = cms.bool(False),
        eid_input_name = cms.string('input'),
        eid_min_cluster_energy = cms.double(1),
        eid_n_clusters = cms.int32(10),
        eid_n_layers = cms.int32(50),
        eid_output_name_energy = cms.string('output/regressed_energy'),
        eid_output_name_id = cms.string('output/id_probabilities'),
        kernelDensityFactor = cms.vdouble(0.1, 0.307354051, 0.4),
        minNumLayerCluster = cms.vint32(4, 4, 3),
        nearestHigherOnSameLayer = cms.bool(False),
        outlierMultiplier = cms.vdouble(1.0, 1.0, 1.37603455),
        rescaleDensityByZ = cms.bool(False),
        type = cms.string('CLUE3D'),
        useAbsoluteProjectiveScale = cms.bool(True),
        useClusterDimensionXY = cms.bool(False)
    )


    process.ticlIterationsTask = cms.Task(
        ticlCLUE3DEMStepTask,
    )


    process.mtdSoA = _mtdSoAProducer.clone()
    process.mtdSoATask = cms.Task(process.mtdSoA)

    #process.ticlTracksterLinks = _tracksterLinksProducer.clone()
    process.ticlTracksterLinks = cms.EDProducer("TracksterLinksProducer",
    detector = cms.string('HGCAL'),
    layer_clusters = cms.InputTag("hgcalMergeLayerClusters"),
    layer_clustersTime = cms.InputTag("hgcalMergeLayerClusters","timeLayerCluster"),
    linkingPSet = cms.PSet(
        algo_verbosity = cms.int32(0),
        alignement_projective_th = cms.double(9.96784068),
        dot_prod_th = cms.double(0.912697019),
        max_distance_closest_points = cms.double(34.7724337),
        max_z_distance_closest_ponts = cms.double(15.6415521),
        min_distance_z = cms.double(19.5239939),
        min_num_lcs = cms.uint32(4),
        min_trackster_energy = cms.double(8.46463942),
        pca_quality_th = cms.double(0.909878691),
        track_time_quality_threshold = cms.double(0.5),
        type = cms.string('Skeletons'),
        wind = cms.double(0.753263398)
    ),
    mightGet = cms.optional.untracked.vstring,
    original_masks = cms.VInputTag("hgcalMergeLayerClusters:InitialLayerClustersMask"),
    propagator = cms.string('PropagatorWithMaterial'),
    tracksters_collections = cms.VInputTag("ticlTrackstersCLUE3DEM")
)
    process.ticlTracksterLinksTask = cms.Task(process.ticlTracksterLinks)
    

    process.ticlCandidate = _ticlCandidateProducer.clone()
    process.ticlCandidateTask = cms.Task(process.ticlCandidate)

    process.tracksterSimTracksterAssociationLinkingbyCLUE3DEM = _tracksterSimTracksterAssociationLinkingbyCLUE3D.clone(
        label_tst = cms.InputTag("ticlTrackstersCLUE3DEM")
        )
    process.tracksterSimTracksterAssociationPRbyCLUE3DEM = _tracksterSimTracksterAssociationPRbyCLUE3D.clone(
        label_tst = cms.InputTag("ticlTrackstersCLUE3DEM")
        )
    #process.tracksterSimTracksterAssociationLinkingbyCLUE3DHAD = _tracksterSimTracksterAssociationLinkingbyCLUE3D.clone(
    #    label_tst = cms.InputTag("ticlTrackstersCLUE3DHAD")
    #    )
    #process.tracksterSimTracksterAssociationPRbyCLUE3DHAD = _tracksterSimTracksterAssociationPRbyCLUE3D.clone(
    #    label_tst = cms.InputTag("ticlTrackstersCLUE3DHAD")
    #    )

    #process.mergedTrackstersProducer = _mergedTrackstersProducer.clone()    

    process.tracksterSimTracksterAssociationLinkingbyCLUE3D = _tracksterSimTracksterAssociationLinkingbyCLUE3D.clone(
        label_tst = cms.InputTag("ticlTrackstersCLUE3DEM")
        )
    process.tracksterSimTracksterAssociationPRbyCLUE3D = _tracksterSimTracksterAssociationPRbyCLUE3D.clone(
        label_tst = cms.InputTag("ticlTrackstersCLUE3DEM")
        )
    process.iterTICLTask = cms.Task(process.ticlLayerTileTask,
                                     process.mtdSoATask,
                                     process.ticlIterationsTask,
                                     process.ticlTracksterLinksTask,
                                     process.ticlCandidateTask)
    process.particleFlowClusterHGCal.initialClusteringStep.tracksterSrc = "ticlCandidate"
    process.globalrecoTask.remove(process.ticlTrackstersMerge)

    process.tracksterSimTracksterAssociationLinking.label_tst = cms.InputTag("ticlCandidate")
    process.tracksterSimTracksterAssociationPR.label_tst = cms.InputTag("ticlCandidate")

    process.tracksterSimTracksterAssociationLinkingPU.label_tst = cms.InputTag("ticlTracksterLinks")
    process.tracksterSimTracksterAssociationPRPU.label_tst = cms.InputTag("ticlTracksterLinks")
    process.mergeTICLTask = cms.Task()
    process.pfTICL = _pfTICLProducerV5.clone()
    process.hgcalAssociators = cms.Task(process.lcAssocByEnergyScoreProducer, process.layerClusterCaloParticleAssociationProducer,
                            process.scAssocByEnergyScoreProducer, process.layerClusterSimClusterAssociationProducer,
                            process.lcSimTSAssocByEnergyScoreProducer, process.layerClusterSimTracksterAssociationProducer,
                            process.simTsAssocByEnergyScoreProducer,  process.simTracksterHitLCAssociatorByEnergyScoreProducer,
                            process.tracksterSimTracksterAssociationLinking, process.tracksterSimTracksterAssociationPR,
                            process.tracksterSimTracksterAssociationLinkingbyCLUE3D, process.tracksterSimTracksterAssociationPRbyCLUE3D,
                            process.tracksterSimTracksterAssociationLinkingbyCLUE3DEM, process.tracksterSimTracksterAssociationPRbyCLUE3DEM,
                            #process.tracksterSimTracksterAssociationLinkingbyCLUE3DHAD, process.tracksterSimTracksterAssociationPRbyCLUE3DHAD,
                            process.tracksterSimTracksterAssociationLinkingPU, process.tracksterSimTracksterAssociationPRPU
                            )

    process.hgcalValidatorv5 = hgcalValidatorv5.clone(
        ticlTrackstersMerge = cms.InputTag("ticlCandidate"),
        trackstersclue3d = cms.InputTag("ticlTrackstersCLUE3DEM")
    )
    process.hgcalValidatorSequence = cms.Sequence(process.hgcalValidatorv5)
    process.hgcalValidation = cms.Sequence(process.hgcalSimHitValidationEE+process.hgcalSimHitValidationHEF+process.hgcalSimHitValidationHEB+process.hgcalDigiValidationEE+process.hgcalDigiValidationHEF+process.hgcalDigiValidationHEB+process.hgcalRecHitValidationEE+process.hgcalRecHitValidationHEF+process.hgcalRecHitValidationHEB+process.hgcalHitValidationSequence+process.hgcalValidatorSequence+process.hgcalTiclPFValidation+process.hgcalPFJetValidation)
    process.globalValidationHGCal = cms.Sequence(process.hgcalValidation)
    process.validation_step9 = cms.EndPath(process.globalValidationHGCal)
    if(enableDumper):
        process.ticlDumper = ticlDumper.clone(
            saveLCs=True,
            saveCLUE3DTracksters=True,
            saveTrackstersMerged=True,
            saveSimTrackstersSC=True,
            saveSimTrackstersCP=True,
            saveTICLCandidate=True,
            saveSimTICLCandidate=True,
            saveTracks=True,
            saveAssociations=True,
            trackstersclue3d = cms.InputTag('ticlTrackstersCLUE3DEM'), 
            ticlcandidates = cms.InputTag("ticlCandidate"),
            trackstersmerged = cms.InputTag("ticlCandidate"),
            trackstersInCand = cms.InputTag("ticlCandidate")
        )
        process.TFileService = cms.Service("TFileService",
                                           fileName=cms.string("histo.root")
                                           )
        process.FEVTDEBUGHLToutput_step = cms.EndPath(
            process.FEVTDEBUGHLToutput + process.ticlDumper)


    process = customiseForTICLv5EventContent(process)

    return process
