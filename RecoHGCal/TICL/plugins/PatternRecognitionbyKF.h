// Author: Mark Matthewman - mark.matthewman@cern.ch
// Date: 04/2021

#ifndef __RecoHGCal_TICL_PRbyKF_H__
#define __RecoHGCal_TICL_PRbyKF_H__
#include <memory>  // unique_ptr

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"

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


namespace ticl {
  template <typename TILES>
  class PatternRecognitionbyKF final : public PatternRecognitionAlgoBaseT<TILES> {
  public:
    PatternRecognitionbyKF(const edm::ParameterSet& conf, edm::ConsumesCollector);
    ~PatternRecognitionbyKF() override = default;

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
    int rescaleFTSError_;
    uint32_t geomCacheId_;


    // Instance Variables
    hgcal::RecHitTools rhtools_;
    std::map<DetId, const HGCRecHit*> hitMap;
    std::map<DetId,LocalError> lerr;
    std::vector<HGCDiskGeomDet *> disksPos_, disksNeg_;

    float etaBinSize = (TILES::constants_type_t::maxEta - TILES::constants_type_t::minEta)/TILES::constants_type_t::nEtaBins;
    float phiBinSize = 2*M_PI/TILES::constants_type_t::nPhiBins;
    int nPhiBin = TILES::constants_type_t::nPhiBins;
    int nEtaBin = TILES::constants_type_t::nEtaBins;

    //Member Functions
    void dumpTiles(const TILES&) const;
    void makeDisks(int subdet, const CaloGeometry* geom_);
    void addDisk(HGCDiskGeomDet *disk) { 
      (disk->zside() > 0 ? disksPos_ : disksNeg_).push_back(disk);
    }
    std::vector<TrajectoryMeasurement> measurements(const TrajectoryStateOnSurface &tsos, 
      const MeasurementEstimator &mest, 
      const TILES &tiles, 
      int depth);
    template<class Start>
    std::vector<TempTrajectory> advanceOneLayer(const Start &start, 
      const HGCDiskGeomDet * disk, 
      const std::vector<HGCDiskGeomDet *>  &disks, 
      const TILES &tiles,
      PropagationDirection direction, 
      bool &isSilicon,
      TempTrajectory traj);
    const HGCDiskGeomDet * nextDisk(const HGCDiskGeomDet * from, 
      PropagationDirection direction, 
      const std::vector<HGCDiskGeomDet *> &vec,
      bool isSilicon) const;
    const HGCDiskGeomDet * switchDisk(const HGCDiskGeomDet * from, 
      const std::vector<HGCDiskGeomDet *> &vec,
      bool isSilicon) const;
    virtual void fillHitMap(std::map<DetId, const HGCRecHit*>& hitMap,
        const HGCRecHitCollection& recHitsEE,
        const HGCRecHitCollection& recHitsFH,
        const HGCRecHitCollection& recHitsBH) const;
    void calculateLocalError(DetId id,
        const HGCalGeometry* hgcalgeom);
    void init(const edm::Event& evt, const edm::EventSetup& es);
  };
}  // namespace ticl
#endif
