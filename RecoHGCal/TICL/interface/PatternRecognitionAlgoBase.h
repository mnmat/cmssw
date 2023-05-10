// Author: Felice Pantaleo - felice.pantaleo@cern.ch
// Date: 09/2018

#ifndef RecoHGCal_TICL_PatternRecognitionAlgoBase_H__
#define RecoHGCal_TICL_PatternRecognitionAlgoBase_H__

#include <memory>
#include <vector>
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/HGCalReco/interface/KFHit.h"
#include "DataFormats/HGCalReco/interface/TICLLayerTile.h"
#include "DataFormats/HGCalReco/interface/TICLSeedingRegion.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "RecoHGCal/TICL/interface/GlobalCache.h"
#include "RecoHGCal/TICL/interface/commons.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"

namespace edm {
  class Event;
  class EventSetup;
}  // namespace edm

namespace ticl {
  template <typename TILES>
  class PatternRecognitionAlgoBaseT {
  public:
    PatternRecognitionAlgoBaseT(const edm::ParameterSet& conf, edm::ConsumesCollector)
        : algo_verbosity_(conf.getParameter<int>("algo_verbosity")) {}
    virtual ~PatternRecognitionAlgoBaseT(){};

    struct Inputs {
      const edm::Event& ev;
      const edm::EventSetup& es;
      const std::vector<reco::CaloCluster>& layerClusters;
      const std::vector<float>& mask;
      const edm::ValueMap<std::pair<float, float>>& layerClustersTime;
      const TILES& tiles;
      const std::vector<TICLSeedingRegion>& regions;
      const tensorflow::Session* tfSession;

      Inputs(const edm::Event& eV,
             const edm::EventSetup& eS,
             const std::vector<reco::CaloCluster>& lC,
             const std::vector<float>& mS,
             const edm::ValueMap<std::pair<float, float>>& lT,
             const TILES& tL,
             const std::vector<TICLSeedingRegion>& rG,
             const tensorflow::Session* tS)
          : ev(eV), es(eS), layerClusters(lC), mask(mS), layerClustersTime(lT), tiles(tL), regions(rG), tfSession(tS) {}
    };
    
    // (mmatthew): makeTrajectories used only by PatternRecognitionByKF. 
    // TODO: combine with makeTracksters as pure virtual function
    virtual void makeTrajectories(const Inputs& input,
                                std::vector<KFHit>& kfhits,
                                std::vector<KFHit>& prophits,
                                float& abs_fail){};

    virtual void makeTracksters(const Inputs& input,
                                std::vector<Trackster>& result,
                                std::unordered_map<int, std::vector<int>>& seedToTracksterAssociation){};

  protected:
    int algo_verbosity_;
  };
}  // namespace ticl

#endif
