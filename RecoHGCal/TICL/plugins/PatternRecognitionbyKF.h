// Author: Marco Rovere - marco.rovere@cern.ch
// Date: 04/2021

#ifndef __RecoHGCal_TICL_PRbyKF_H__
#define __RecoHGCal_TICL_PRbyKF_H__
#include <memory>  // unique_ptr
#include "RecoHGCal/TICL/interface/PatternRecognitionAlgoBase.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/PatternTools/interface/TempTrajectory.h"

#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/PtrVector.h"
#include "DataFormats/Common/interface/RefProd.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefVector.h"

#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include "Geometry/CommonTopologies/interface/HGCDiskGeomDet.h"
#include "DataFormats/HGCTrackingRecHit/interface/HGCTrackingRecHit.h"

#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimatorBase.h"
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h"

#include "TrackingTools/PatternTools/interface/TrajectoryStateUpdator.h"

#include "Geometry/HGCalCommonData/interface/HGCalDDDConstants.h"

namespace ticl {
  template <typename TILES>
  class PatternRecognitionbyKF final : public PatternRecognitionAlgoBaseT<TILES> {
  public:
    PatternRecognitionbyKF(const edm::ParameterSet& conf, edm::ConsumesCollector);
    ~PatternRecognitionbyKF() override = default;

    template<class Start>
    std::vector<TempTrajectory>
    advanceOneLayer(const Start &start, 
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
        const HGCalDDDConstants* ddd);

    void makeTracksters(const typename PatternRecognitionAlgoBaseT<TILES>::Inputs& input,
                        std::vector<Trackster>& result,
                        std::unordered_map<int, std::vector<int>>& seedToTracksterAssociation) override;


    void makeTracksters_verbose(const typename PatternRecognitionAlgoBaseT<TILES>::Inputs& input,
                        std::vector<Trackster>& result,
                        std::vector<GlobalPoint>& points_kf,
                        std::vector<GlobalPoint>& points_prop,
                        std::unordered_map<int, std::vector<int>>& seedToTracksterAssociation) override;


/*

    void makeTracksters(const typename PatternRecognitionAlgoBaseKFT<TILES>::Inputs& input,
                    std::vector<Trackster>& result,
                    std::vector<GlobalPoint>& points,
                    std::unordered_map<int, std::vector<int>>& seedToTracksterAssociation);

*/
    void energyRegressionAndID(const std::vector<reco::CaloCluster>& layerClusters,
                               const tensorflow::Session*,
                               std::vector<Trackster>& result);

    static void fillPSetDescription(edm::ParameterSetDescription& iDesc);

  private:

    // Declarations for Constructor


    void dumpTiles(const TILES&) const;

    edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeomToken_;
    const std::string propName_;
    const std::string propNameOppo_;
    //const std::string propNameRK_;
    edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> bfieldtoken_;
    edm::ESGetToken<Propagator, TrackingComponentsRecord> propagatortoken_;
    edm::ESGetToken<Propagator, TrackingComponentsRecord> propagatorOppoToken_;
    edm::ESGetToken<Chi2MeasurementEstimatorBase, TrackingComponentsRecord> estimatorToken_;
    edm::ESGetToken<TrajectoryStateUpdator, TrackingComponentsRecord> updatorToken_;
    // edm::ESGetToken<Propagator, TrackingComponentsRecord> propagatortokenRK_;
    edm::EDGetTokenT<reco::TrackCollection> trackToken_;
    edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsEEToken_; 
    edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsFHToken_;
    edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsBHToken_;

    const std::string eidInputName_;
    const std::string eidOutputNameEnergy_;
    const std::string eidOutputNameId_;
    const float eidMinClusterEnergy_;
    const int eidNLayers_;
    const int eidNClusters_;
    const std::string materialbudget_;

    edm::ESHandle<MagneticField> bfield_;
    edm::ESHandle<Propagator> propagator_;
    edm::ESHandle<Chi2MeasurementEstimatorBase> estimator_;
    edm::ESHandle<TrajectoryStateUpdator> updator_;
    edm::ESHandle<Propagator> propagatorOppo_;

    //edm::ESHandle<Propagator> propagatorRK_;

    std::map<std::string, float> xi_;

    hgcal::RecHitTools rhtools_;
    tensorflow::Session* eidSession_;

    static const int eidNFeatures_ = 3;
    std::map<DetId, const HGCRecHit*> hitMap;

    //TILE constants

    float etaBinSize = (TILES::constants_type_t::maxEta - TILES::constants_type_t::minEta)/TILES::constants_type_t::nEtaBins;
    float phiBinSize = 2*M_PI/TILES::constants_type_t::nPhiBins;
    int nPhiBin = TILES::constants_type_t::nPhiBins;
    int nEtaBin = TILES::constants_type_t::nEtaBins;

    std::map<DetId,std::pair<float,float>> lerr;

    std::vector<HGCDiskGeomDet *> disksPos_, disksNeg_;

    void makeDisks(int subdet, int disks, const CaloGeometry* geom_);
    void addDisk(HGCDiskGeomDet *disk) { 
      (disk->zside() > 0 ? disksPos_ : disksNeg_).push_back(disk);
    }

    std::vector<TrajectoryMeasurement> measurements(const TrajectoryStateOnSurface &tsos, const MeasurementEstimator &mest, const TILES &tiles, int depth);


  };

}  // namespace ticl
#endif
