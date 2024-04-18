import FWCore.ParameterSet.Config as cms

from RecoHGCal.TICL.TICLSeedingRegions_cff import ticlSeedingTrk
from RecoHGCal.TICL.trackstersProducer_cfi import trackstersProducer as _trackstersProducer
from RecoHGCal.TICL.filteredLayerClustersProducer_cfi import filteredLayerClustersProducer as _filteredLayerClustersProducer

from RecoHGCal.TICL.multiClustersFromTrackstersProducer_cfi import multiClustersFromTrackstersProducer as _multiClustersFromTrackstersProducer
from RecoHGCal.TICL.ticlLayerTileProducer_cfi import ticlLayerTileProducer

# CLUSTER FILTERING/MASKING
# Layer Clusters not used by Kalman Filter implementation but necessary for TracksterProducer. Masking strategy a work in progress!

filteredLayerClustersKalmanFilter = _filteredLayerClustersProducer.clone(
    clusterFilter = "ClusterFilterByAlgoAndSize",
    min_cluster_size = 2, # inclusive
    algo_number = [8],
    iteration_label = "KalmanFilter",
)

# PATTERN RECOGNITION

ticlTrackstersKalmanFilter = _trackstersProducer.clone(
    filtered_mask = "filteredLayerClustersKalmanFilter:KalmanFilter",
    layer_clusters_tiles = "ticlRecHitTile",
    seeding_regions = "ticlSeedingTrk",
    itername = "KalmanFilter",
    patternRecognitionBy = "KalmanFilter",
    pluginPatternRecognitionByKalmanFilter = dict (
        rescaleFTSError = 1., # used to rescale the Error of the last FTS of the Tracker which is propagated to the first layer of HGCAL
        scaleWindow = 1.,
    )
)

HGCTrackerESProducer = cms.ESProducer("HGCTrackerESProducer")

ticlKalmanFilterStepTask = cms.Task(ticlSeedingTrk
    ,filteredLayerClustersKalmanFilter
    ,HGCTrackerESProducer
    ,ticlTrackstersKalmanFilter)
