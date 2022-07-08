#ifndef RecoHGCal_HGCTracking_HGCDiskGeomDet_h
#define RecoHGCal_HGCTracking_HGCDiskGeomDet_h

/// Class corresponding to one layer of HGC

#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCScintillatorDetId.h"
#include "FWCore/Utilities/interface/Exception.h"

class HGCDiskGeomDet : public GeomDet {
    public:
        HGCDiskGeomDet(int subdet, int zside, int layer, float z, float rmin, float rmax, float radlen, float xi) ;

        int subdet() const { return subdet_; }
        int zside() const { return zside_; }
        int layer() const { return layer_; }
        float rmin() const { return rmin_; }
        float rmax() const { return rmax_; }
        bool isSilicon() const {
            if(subdet_ == 8 || subdet_ == 9){
                return true;
            }
            else if (subdet_ == 10){
                return false;
            }
            else throw cms::Exception("LogicError", "Subdetector not defined");
        }

    protected:
        const int subdet_, zside_, layer_;
        const float rmin_, rmax_;
};

#endif