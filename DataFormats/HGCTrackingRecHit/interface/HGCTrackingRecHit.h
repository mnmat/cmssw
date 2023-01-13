#ifndef RecoHGCal_HGCTracking_HGCTrackingRecHit_h
#define RecoHGCal_HGCTracking_HGCTrackingRecHit_h

/// Basic template class for a RecHit wrapping a Ref to an object

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