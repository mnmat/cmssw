// Author: Mark Matthewman - mark.matthewman@cern.ch
// Date: 07/2026

#ifndef __RecoHGCal_TICL_HGCTrackingbyKalmanFilter_H__
#define __RecoHGCal_TICL_HGCTrackingbyKalmanFilter_H__
#include <memory>  // unique_ptr

#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/HGCalReco/interface/HGCTrackingRecHit.h"

#include "Geometry/CommonTopologies/interface/HGCDiskGeomDet.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"


#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "RecoHGCal/TICL/interface/HGCTrackingAlgoBase.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/PatternTools/interface/TrajectoryStateUpdator.h"
#include "TrackingTools/PatternTools/interface/TempTrajectory.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimatorBase.h"
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h"

#include "HGCTracker.h"

namespace ticl {
  template <typename TILES>
  class HGCTrackingbyKalmanFilter final : public HGCTrackingAlgoBaseT<TILES> {
  public:
    HGCTrackingbyKalmanFilter(const edm::ParameterSet& conf, edm::ConsumesCollector);
    ~HGCTrackingbyKalmanFilter() override = default;

    void makeTrajectories(const typename HGCTrackingAlgoBaseT<TILES>::Inputs& input,
                        std::vector<KFHit>& kfhits,
                        std::vector<reco::Track>& tracks,
                        std::vector<reco::TrackExtra>& trackExtras,
                        TrackingRecHitCollection& trackingRecHitCollection) override;

    static void fillPSetDescription(edm::ParameterSetDescription& iDesc);
    void setGeometry(hgcal::RecHitTools const& rhtools) override {};

  private:
    edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeomToken_;
    const std::string propName_;
    const std::string propNameOppo_;
    edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> bfieldtoken_;
    edm::ESGetToken<Propagator, TrackingComponentsRecord> propagatortoken_;
    edm::ESGetToken<Propagator, TrackingComponentsRecord> propagatorOppoToken_;
    edm::ESGetToken<Chi2MeasurementEstimatorBase, TrackingComponentsRecord> estimatorToken_;
    edm::ESGetToken<TrajectoryStateUpdator, TrackingComponentsRecord> updatorToken_;
    edm::EDGetTokenT<reco::TrackCollection> trackToken_;
    edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsEEToken_; 
    edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsFHToken_;
    edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsBHToken_;

    edm::ESHandle<MagneticField> bfield_;
    edm::ESHandle<Propagator> propagator_;
    edm::ESHandle<Chi2MeasurementEstimatorBase> estimator_;
    edm::ESHandle<TrajectoryStateUpdator> updator_;
    edm::ESHandle<Propagator> propagatorOppo_;
    edm::ESGetToken<HGCDiskGeomDetVector,CaloGeometryRecord> diskToken_;
    edm::ESGetToken<HGCTracker,CaloGeometryRecord> hgcTrackerToken_;
    edm::Handle<HGCRecHitCollection> ee_hits;
    edm::Handle<HGCRecHitCollection> fh_hits;
    edm::Handle<HGCRecHitCollection> bh_hits;
    
    double rescaleFTSError_;
    double scaleWindow_;

    bool standalonePropagator_;
    bool doBackwardPropagation_;
    uint64_t geomCacheId_;
    int trackId;
    int evtId;

    // Instance Variables
    hgcal::RecHitTools rhtools_;
    std::vector<std::pair<const HGCRecHit*, int>> recHitCollection;
    const HGCTracker* hgcTracker_;

    enum TColl{
      HGCEERecHits,
      HGCHEFRecHits,
      HGCHEBRecHits,    
    };

    static constexpr float etaBinSize = (TILES::constants_type_t::maxEta - TILES::constants_type_t::minEta)/TILES::constants_type_t::nEtaBins;
    static constexpr float phiBinSize = 2*M_PI/TILES::constants_type_t::nPhiBins;
    static constexpr int nPhiBin = TILES::constants_type_t::nPhiBins;
    static constexpr int nEtaBin = TILES::constants_type_t::nEtaBins;

    //Member Functions
    std::pair<float,float> covarianceTransform(const TrajectoryStateOnSurface &tsos);
    void dumpTiles(const TILES&) const;
    std::vector<std::shared_ptr<HGCTrackingRecHit>> measurements(const TrajectoryStateOnSurface &tsos, 
      const TILES &tiles, 
      int depth);
    template<class Start>
    std::vector<TempTrajectory> advanceOneLayer(const Start &start, 
      const HGCDiskLayer * disk,
      const TILES &tiles,
      PropagationDirection direction, 
      bool &isSilicon,
      TempTrajectory traj);
    virtual void mergeRecHitCollections(std::vector<std::pair<const HGCRecHit*, int>>& recHitCollection,
        const HGCRecHitCollection& recHitsEE,
        const HGCRecHitCollection& recHitsFH,
        const HGCRecHitCollection& recHitsBH) const;
    void init(const edm::Event& evt, const edm::EventSetup& es);
  };
}  // namespace ticl

#endif