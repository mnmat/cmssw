#ifndef RecoHGCal_HGCTracking_HGCTrackingRecHit_h
#define RecoHGCal_HGCTracking_HGCTrackingRecHit_h

/// Wrapper for TrackingRecHits for PatternRecognitionByKF

#include <cassert>
#include "DataFormats/TrackingRecHit/interface/RecHit2DLocalPos.h"

class HGCTrackingRecHit : public RecHit2DLocalPos {
    public:
        HGCTrackingRecHit() {}
        HGCTrackingRecHit(DetId id, const LocalPoint &pos, const LocalError &err):
            RecHit2DLocalPos(id),
            pos_(pos), err_(err) {}

        virtual HGCTrackingRecHit * clone() const override { return new HGCTrackingRecHit(*this); }
        virtual LocalPoint localPosition() const override { return pos_; }
        virtual LocalError localPositionError() const override { return err_; }

    protected:
        LocalPoint pos_;
        LocalError err_;
        float energy_;
};

#endif