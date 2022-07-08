#ifndef RecoHGCal_HGCTracking_HGCTrackingRecHit_h
#define RecoHGCal_HGCTracking_HGCTrackingRecHit_h

/// Wrapper for TrackingRecHits for PatternRecognitionByKF

#include <cassert>
#include "DataFormats/TrackingRecHit/interface/RecHit2DLocalPos.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"


class HGCTrackingRecHit : public RecHit2DLocalPos {
    public:
        HGCTrackingRecHit() {}
        HGCTrackingRecHit(DetId id, const HGCRecHitRef &objref, const LocalPoint &pos, const LocalError &err):
            RecHit2DLocalPos(id),
            objref_(objref),
            pos_(pos), err_(err) {}

        virtual HGCTrackingRecHit* clone() const override { return new HGCTrackingRecHit(*this); }
        virtual LocalPoint localPosition() const override { return pos_; }
        virtual LocalError localPositionError() const override { return err_; }

        const HGCRecHitRef & objRef() const { return objref_; }

        bool sharesInput( const TrackingRecHit* other, SharedInputType what) const {
            if (typeid(*other) == typeid(HGCTrackingRecHit)) {
                return objRef() == (static_cast<const HGCTrackingRecHit *>(other))->objRef();
            } else {
                return false;
            }
        }

    protected:
        HGCRecHitRef objref_;
        LocalPoint pos_;
        LocalError err_;
        float energy_;
};


#endif