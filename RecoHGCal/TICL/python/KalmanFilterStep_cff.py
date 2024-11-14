import FWCore.ParameterSet.Config as cms

from RecoHGCal.TICL.TICLSeedingRegions_cff import ticlSeedingTrk
from RecoHGCal.TICL.trackstersProducer_cfi import trackstersProducer as _trackstersProducer
from RecoHGCal.TICL.filteredLayerClustersProducer_cfi import filteredLayerClustersProducer as _filteredLayerClustersProducer
from RecoHGCal.TICL.ticlLayerTileProducer_cfi import ticlLayerTileProducer
from TrackPropagation.Geant4e.Geant4ePropagator_cfi import Geant4ePropagator

# CLUSTER FILTERING/MASKING
# Layer Clusters not used by Kalman Filter implementation but necessary for TracksterProducer. Masking strategy a work in progress!

filteredLayerClustersKalmanFilter = _filteredLayerClustersProducer.clone(
    clusterFilter = "ClusterFilterByAlgoAndSize",
    min_cluster_size = 2, # inclusive
    algo_number = [8],
    iteration_label = "KalmanFilter",
)


# PATTERN RECOGNITION

from SimG4Core.Application.g4SimHits_cfi import g4SimHits as _g4SimHits

geoprotest = cms.EDProducer("GeometryProducer",
     GeoFromDD4hep = cms.bool(False),
     UseMagneticField = cms.bool(True),
     UseSensitiveDetectors = cms.bool(False),
     MagneticField =  _g4SimHits.MagneticField.clone()
)

#from Configuration.ProcessModifiers.dd4hep_cff import dd4hep
#dd4hep.toModify(geoprotest, GeoFromDD4hep = True )

Geant4ePropagatorTest = Geant4ePropagator.clone(
    ComponentName=cms.string("Geant4ePropagatorTest")
)

ticlTrackstersKalmanFilter = _trackstersProducer.clone(
    filtered_mask = "filteredLayerClustersKalmanFilter:KalmanFilter",
    layer_clusters_tiles = "ticlRecHitTile",
    seeding_regions = "ticlSeedingTrk",
    itername = "KalmanFilter",
    patternRecognitionBy = "KalmanFilter",
    pluginPatternRecognitionByKalmanFilter = dict (
        rescaleFTSError = 2., # used to rescale the Error of the last FTS of the Tracker which is propagated to the first layer of HGCAL
        scaleWindow = 1.,
        #propagator = "RungeKuttaTrackerPropagator",
        propagator = "Geant4ePropagatorTest",
    )
)

ticlTrackstersStandalonePropagator = _trackstersProducer.clone(
    filtered_mask = "filteredLayerClustersKalmanFilter:KalmanFilter",
    layer_clusters_tiles = "ticlRecHitTile",
    seeding_regions = "ticlSeedingTrk",
    itername = "KalmanFilter",
    patternRecognitionBy = "KalmanFilter",
    pluginPatternRecognitionByKalmanFilter = dict (
        rescaleFTSError = 2., # used to rescale the Error of the last FTS of the Tracker which is propagated to the first layer of HGCAL
        scaleWindow = 1.,
        standalonePropagator = True,
        #propagator = "RungeKuttaTrackerPropagator",
        propagator = "Geant4ePropagatorTest",
    )
)

HGCTrackerESProducer = cms.ESProducer("HGCTrackerESProducer",
    radlen = cms.vdouble(1.53770787, 0.71064359, 1.45345887, 0.57315113, 1.02882455,
        0.92384098, 1.15461784, 0.72404336, 0.9948446 , 1.0107427 ,
        1.05947235, 0.91730167, 1.20028302, 0.6703572 , 0.98144224,
        1.01024202, 1.08452792, 0.86282149, 1.53923452, 0.99185102,
        1.67874405, 0.70709974, 1.63824099, 0.97162878, 1.74571227,
        0.69011827, 2.92834302, 3.01147101, 3.0583451 , 3.12601533,
        2.85205937, 2.95217992, 3.14263578, 3.07471756, 3.05502943,
        2.82345623, 3.0230636 , 4.29398744, 3.9234094 , 4.27748842,
        3.91229994, 4.23728221, 4.02845205, 4.21537293, 4.32452121,
        3.83363941, 4.32332509),
    xi = cms.vdouble(
        0.00264665, 0.00050171, 0.00081145, 0.0003883 , 0.00049233,
        0.00066116, 0.00059059, 0.00050953, 0.0004874 , 0.00069975,
        0.00051153, 0.00065396, 0.00062511, 0.00046753, 0.00048085,
        0.00070629, 0.00052536, 0.00061608, 0.00068004, 0.0006793 ,
        0.00079585, 0.0005112 , 0.00073527, 0.00067903, 0.00081588,
        0.00049695, 0.00292419, 0.0029901 , 0.00304991, 0.00309913,
        0.00283049, 0.00295079, 0.00310096, 0.00307401, 0.00304154,
        0.00278516, 0.00302603, 0.00425792, 0.00389764, 0.00425274,
        0.0038517 , 0.00426734, 0.00399329, 0.00420363, 0.00429176,
        0.00378711, 0.004295657)
    )

ticlKalmanFilterStepTask = cms.Task(
    geoprotest
    ,ticlSeedingTrk
    ,filteredLayerClustersKalmanFilter
    ,HGCTrackerESProducer
    ,ticlTrackstersKalmanFilter
    ,ticlTrackstersStandalonePropagator)
