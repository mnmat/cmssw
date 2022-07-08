// Author: Mark Matthewman - mark.matthewman@cern.ch
// Date: 04/2023

// Note: KFHit is a temporary struct provided to output the TSOS produced by the PatternRecognitionbyKalmanFilter in HGCAL.
// Future iterations of the algorithm won't produce KFHits and then this struct is obsolete.

#ifndef DataFormats_HGCalReco_KFHit_h
#define DataFormats_HGCalReco_KFHit_h

#include <array>
#include <vector>
#include "DataFormats/Provenance/interface/ProductID.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include <Eigen/Core>

struct KFHit
{  
    explicit KFHit(TrajectoryStateOnSurface tsos, uint32_t detid, reco::Track track, int trackId, int layer): 
        center(tsos.globalPosition()),
        localError(tsos.localError().matrix()),
        cartesianError(tsos.cartesianError().matrix()),
        curvilinearError(tsos.curvilinearError().matrix()),
        xx(tsos.localError().positionError().xx()),
        xy(tsos.localError().positionError().xy()),
        yy(tsos.localError().positionError().yy()),
        charge(tsos.charge()),
        detid(detid),
        eta(tsos.localMomentum().eta()),
        theta(tsos.localMomentum().theta()),
        trackId(trackId),
        trackCharge(track.charge()),
        trackMomentum(track.p()),
        trackQuality(track.qualityMask()),
        trackChi2(track.chi2()),
        trackValidFraction(track.validFraction()),
        trackQOverP(track.qoverp()),
        layer(layer)
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
    float eta;
    float theta;
    int trackId;
    int trackCharge;
    float trackMomentum;
    float trackQuality;
    float trackChi2;
    float trackValidFraction;
    float trackQOverP;
    int layer;
};
#endif
