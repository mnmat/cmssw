#include "RecoHGCal/TICL/plugins/HGCTrackingPluginFactory.h"
#include "HGCTrackingbyKalmanFilter.h"
#include "FWCore/ParameterSet/interface/ValidatedPluginFactoryMacros.h"
#include "FWCore/ParameterSet/interface/ValidatedPluginMacros.h"
EDM_REGISTER_VALIDATED_PLUGINFACTORY(HGCTrackingFactory, "HGCTrackingFactory");

DEFINE_EDM_VALIDATED_PLUGIN(HGCTrackingFactory, ticl::HGCTrackingbyKalmanFilter<TICLLayerTiles>, "KalmanFilter");

