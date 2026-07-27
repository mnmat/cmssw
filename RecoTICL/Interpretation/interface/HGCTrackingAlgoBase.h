#ifndef RecoHGCal_TICL_HGCTrackingAlgoBase_h
#define RecoHGCal_TICL_HGCTrackingAlgoBase_h

#include <memory>
#include <vector>
#include <functional>
#include <algorithm>

#include "DataFormats/HGCalReco/interface/KFHit.h"
#include "DataFormats/HGCalReco/interface/TICLLayerTile.h"
#include "DataFormats/HGCalReco/interface/TICLSeedingRegion.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

namespace edm {
  class Event;
  class EventSetup;
}  // namespace edm

namespace ticl{
  template <typename TILES>
  class HGCTrackingAlgoBaseT {
  public:
    HGCTrackingAlgoBaseT(const edm::ParameterSet& conf, edm::ConsumesCollector)
        : algo_verbosity_(conf.getParameter<int>("algo_verbosity")) {}
    virtual ~HGCTrackingAlgoBaseT() {};

    struct Inputs {
      edm::Event& ev;
      const edm::EventSetup& es;
      const TILES& tiles;
      const std::vector<TICLSeedingRegion>& regions;
      Inputs(edm::Event& eV,
             const edm::EventSetup& eS,
             const TILES& tL,
             const std::vector<TICLSeedingRegion>& rG)
          : ev(eV), es(eS), tiles(tL), regions(rG) {}
    };

    virtual void makeTrajectories(const Inputs& input,
                                std::vector<KFHit>& kfhits,
                                std::vector<reco::Track>& tracks,
                                std::vector<reco::TrackExtra>& trackExtras,
                                TrackingRecHitCollection& trackingRecHitCollection) = 0;

    virtual void setGeometry(hgcal::RecHitTools const& rhtools) = 0;

  protected:
    int algo_verbosity_;

    hgcal::RecHitTools const* rhtools_ = nullptr; // non-owning, set in the beginRun()
  };
} // namespace ticl

#endif