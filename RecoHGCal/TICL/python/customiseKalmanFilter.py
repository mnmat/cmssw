#from RecoHGCal.TICL.iterativeTICL_cff import ticlIterLabels, ticlIterLabelsMerge
from RecoHGCal.TICL.KalmanFilterStep_cff import ticlTrackstersKalmanFilter, ticlTrackstersStandalonePropagator

def customiseKalmanFilter(process,mode):
    process.ticlTrackstersKalmanFilter.pluginPatternRecognitionByKalmanFilter.distanceRequirementMode = mode
    process.ticlTrackstersStandalonePropagator.pluginPatternRecognitionByKalmanFilter.distanceRequirementMode = mode
    process.ticlTrackstersStandalonePropagatorG4e.pluginPatternRecognitionByKalmanFilter.distanceRequirementMode = mode
    process.ticlTrackstersKalmanFilterG4e.pluginPatternRecognitionByKalmanFilter.distanceRequirementMode = mode
    return process

def customiseRescaleFTS(process,rescaleFTS):
    process.ticlTrackstersKalmanFilter.pluginPatternRecognitionByKalmanFilter.rescaleFTSError = rescaleFTS
    process.ticlTrackstersStandalonePropagator.pluginPatternRecognitionByKalmanFilter.rescaleFTSError = rescaleFTS
    process.ticlTrackstersStandalonePropagatorG4e.pluginPatternRecognitionByKalmanFilter.rescaleFTSError = rescaleFTS
    process.ticlTrackstersKalmanFilterG4e.pluginPatternRecognitionByKalmanFilter.rescaleFTSError = rescaleFTS
    return process

def customisePropagator(process,propagator):
    process.ticlTrackstersKalmanFilter.pluginPatternRecognitionByKalmanFilter.propagator = propagator
    process.ticlTrackstersStandalonePropagator.pluginPatternRecognitionByKalmanFilter.propagator = propagator
    return process

