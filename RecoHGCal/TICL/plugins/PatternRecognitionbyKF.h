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

namespace ticl {
  template <typename TILES>
  class PatternRecognitionbyKF final : public PatternRecognitionAlgoBaseT<TILES> {
  public:
    PatternRecognitionbyKF(const edm::ParameterSet& conf, edm::ConsumesCollector);
    ~PatternRecognitionbyKF() override = default;

    const GeomDet * nextDisk(const GeomDet * from, 
                    PropagationDirection direction, 
                    const std::vector<GeomDet *> &vec) const;

    void makeTracksters(const typename PatternRecognitionAlgoBaseT<TILES>::Inputs& input,
                        std::vector<Trackster>& result,
                        std::unordered_map<int, std::vector<int>>& seedToTracksterAssociation) override;


    void makeTracksters_verbose(const typename PatternRecognitionAlgoBaseT<TILES>::Inputs& input,
                        std::vector<Trackster>& result,
                        std::vector<GlobalPoint>& points,
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

    edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeomToken_;
    const std::string propName_;
    edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> bfieldtoken_;
    edm::ESGetToken<Propagator, TrackingComponentsRecord> propagatortoken_;
    edm::EDGetTokenT<reco::TrackCollection> trackToken_;

    const std::string eidInputName_;
    const std::string eidOutputNameEnergy_;
    const std::string eidOutputNameId_;
    const float eidMinClusterEnergy_;
    const int eidNLayers_;
    const int eidNClusters_;

    edm::ESHandle<MagneticField> bfield_;
    edm::ESHandle<Propagator> propagator_;

    std::map<std::string, float> xi_;

    hgcal::RecHitTools rhtools_;
    tensorflow::Session* eidSession_;

    static const int eidNFeatures_ = 3;

    std::vector<GeomDet *> disksPos_, disksNeg_;

    void computeAbsorbers();
    float combinedEdX(float w1, float a, float w2, float b){
      return (w1 * a + w2 * b);
    };
    float combineX0(float w1, float a, float w2, float b){
      float oneOver = (w1 / a + w2 / b);
      return 1./oneOver;
    };
    void makeDisks(int subdet, int disks, const CaloGeometry* geom_);
    void addDisk(GeomDet *disk, int zside){
      (zside > 0 ? disksPos_ : disksNeg_).push_back(disk);
    }


  };

}  // namespace ticl
#endif
