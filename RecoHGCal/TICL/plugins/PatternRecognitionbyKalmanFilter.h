// Author: Mark Matthewman - mark.matthewman@cern.ch
// Date: 04/2021

#ifndef __RecoHGCal_TICL_PRbyKF_H__
#define __RecoHGCal_TICL_PRbyKF_H__
#include <memory>  // unique_ptr

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/HGCalReco/interface/HGCTrackingRecHit.h"

#include "Geometry/CommonTopologies/interface/HGCDiskGeomDet.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"


#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "RecoHGCal/TICL/interface/PatternRecognitionAlgoBase.h"

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
  class PatternRecognitionbyKalmanFilter final : public PatternRecognitionAlgoBaseT<TILES> {
  public:
    PatternRecognitionbyKalmanFilter(const edm::ParameterSet& conf, edm::ConsumesCollector);
    ~PatternRecognitionbyKalmanFilter() override = default;

    void makeTracksters(const typename PatternRecognitionAlgoBaseT<TILES>::Inputs& input,
                        std::vector<Trackster>& result,
                        std::unordered_map<int, std::vector<int>>& seedToTracksterAssociation) override {};

    void makeTrajectories(const typename PatternRecognitionAlgoBaseT<TILES>::Inputs& input,
                        std::vector<KFHit>& kfhits) override;

    static void fillPSetDescription(edm::ParameterSetDescription& iDesc);

  private:
    // Declarations for Constructor
    edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeomToken_;
    const std::vector<double> radlen_;
    const std::vector<double> xi_;
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

    double rescaleFTSError_;
    double scaleWindow_;

    uint32_t geomCacheId_;

    // Instance Variables
    hgcal::RecHitTools rhtools_;
    std::vector<const HGCRecHit*> recHitCollection;
    std::map<DetId,LocalError> lerr;
    const HGCTracker* hgcTracker_;

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
    virtual void mergeRecHitCollections(std::vector<const HGCRecHit*>& recHitCollection,
        const HGCRecHitCollection& recHitsEE,
        const HGCRecHitCollection& recHitsFH,
        const HGCRecHitCollection& recHitsBH) const;
    void init(const edm::Event& evt, const edm::EventSetup& es);
  };
}  // namespace ticl
#endif
