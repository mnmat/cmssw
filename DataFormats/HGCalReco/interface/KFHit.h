// Author: Mark Matthewman - mark.matthewman@cern.ch
// Date: 04/2023

// Note!!!!
// KFHit is a temporary class provided for easier handling of the PatternRecognitionByKF.
// Once the result format for the Algo has been decided, then the KFHit is obsolete and needs to be deleated.

#ifndef DataFormats_HGCalReco_KFHit_h
#define DataFormats_HGCalReco_KFHit_h

#include <array>
#include <vector>
#include "DataFormats/Provenance/interface/ProductID.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include <Eigen/Core>

struct KFHit
{  
    explicit KFHit(TrajectoryStateOnSurface tsos, uint32_t id): 
        center(tsos.globalPosition()),
        localError(tsos.localError().matrix()),
        cartesianError(tsos.cartesianError().matrix()),
        curvilinearError(tsos.curvilinearError().matrix()),
        xx(tsos.localError().positionError().xx()),
        xy(tsos.localError().positionError().xy()),
        yy(tsos.localError().positionError().yy()),
        charge(tsos.charge()),
        detid(id)
        {}
    KFHit(){}

    GlobalPoint center;
    AlgebraicSymMatrix55 localError;
    AlgebraicSymMatrix66 cartesianError;
    AlgebraicSymMatrix55 curvilinearError;
    float xx;
    float xy;
    float yy;
    int charge;
    uint32_t detid;
};
#endif
