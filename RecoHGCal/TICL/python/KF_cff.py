import FWCore.ParameterSet.Config as cms

from RecoHGCal.TICL.TICLSeedingRegions_cff import ticlSeedingGlobal, ticlSeedingGlobalHFNose
from RecoHGCal.TICL.trackstersProducer_cfi import trackstersProducer as _trackstersProducer
from RecoHGCal.TICL.filteredLayerClustersProducer_cfi import filteredLayerClustersProducer as _filteredLayerClustersProducer
from RecoHGCal.TICL.multiClustersFromTrackstersProducer_cfi import multiClustersFromTrackstersProducer as _multiClustersFromTrackstersProducer

# CLUSTER FILTERING/MASKING

filteredLayerClustersKF = _filteredLayerClustersProducer.clone(
    clusterFilter = "ClusterFilterByAlgoAndSize",
    min_cluster_size = 2, # inclusive
    algo_number = 8,
    iteration_label = "KF"
)

# PATTERN RECOGNITION

ticlTrackstersKF = _trackstersProducer.clone(
    filtered_mask = "filteredLayerClustersKF:KF",
    seeding_regions = "ticlSeedingGlobal",
    itername = "KF",
    patternRecognitionBy = "KF",
    pluginPatternRecognitionByKF = dict (
        criticalEtaPhiDistance = 0.025
    )

)

ticlKFTask = cms.Task(ticlSeedingGlobal
    ,filteredLayerClustersKF
    ,ticlTrackstersKF)

