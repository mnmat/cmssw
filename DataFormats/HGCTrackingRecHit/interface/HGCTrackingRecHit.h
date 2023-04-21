#ifndef RecoHGCal_HGCTracking_HGCTrackingRecHit_h
#define RecoHGCal_HGCTracking_HGCTrackingRecHit_h

/// Wrapper for TrackingRecHits for PatternRecognitionByKF

#include <cassert>
#include "DataFormats/TrackingRecHit/interface/RecHit2DLocalPos.h"

class HGCTrackingRecHit : public RecHit2DLocalPos {
    public:
        HGCTrackingRecHit() {}
        HGCTrackingRecHit(DetId id, const LocalPoint &pos, const LocalError &err, const float energy):
            RecHit2DLocalPos(id),
            pos_(pos), err_(err), energy_(energy) {}

        virtual HGCTrackingRecHit * clone() const override { return new HGCTrackingRecHit(*this); }
        virtual LocalPoint localPosition() const override { return pos_; }
        virtual LocalError localPositionError() const override { return err_; }
        float energy() const { return energy_; }

    protected:
        LocalPoint pos_;
        LocalError err_;
        float energy_;
};

#endif