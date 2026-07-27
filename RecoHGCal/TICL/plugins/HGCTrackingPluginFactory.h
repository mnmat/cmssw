#ifndef RecoHGCal_TICL_HGCTrackingPluginFactory_h
#define RecoHGCal_TICL_HGCTrackingPluginFactory_h

#include "FWCore/PluginManager/interface/PluginFactory.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "RecoHGCal/TICL/interface/HGCTrackingAlgoBase.h"

typedef edmplugin::PluginFactory<ticl::HGCTrackingAlgoBaseT<TICLLayerTiles>*(const edm::ParameterSet&,
                                                                                  edm::ConsumesCollector)>
    HGCTrackingFactory;

#endif