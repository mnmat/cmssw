#ifndef RecoHGCal_HGCTracking_HGCDiskGeomDet_h
#define RecoHGCal_HGCTracking_HGCDiskGeomDet_h

// Wrapper class for GeomDet corresponding to one layer of HGCal. Used by PatternRecognitionByKF.

#include "Geometry/CommonTopologies/interface/GeomDet.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/DetId/interface/DetId.h"

class HGCDiskGeomDet : public GeomDet {
    public:
        HGCDiskGeomDet(int subdet, int zside, int layer, float z, float rmin, float rmax, float radlen, float xi) ;

        int subdet() const { return subdet_; }
        int zside() const { return zside_; }
        int layer() const { return layer_; }
        float rmin() const { return rmin_; }
        float rmax() const { return rmax_; }
        bool isSilicon() const { // The enum values were taken from DataFormats/DetId/interface/DetId.h
            if(subdet_ == DetId::HGCalEE || subdet_ == DetId::HGCalHSi){
                return true;
            }
            else if (subdet_ == DetId::HGCalHSc){
                return false;
            }
            else throw cms::Exception("LogicError", "Subdetector not defined");
        }

    protected:
        const int subdet_, zside_, layer_;
        const float rmin_, rmax_;
};


typedef std::pair<std::vector<HGCDiskGeomDet*>,std::vector<HGCDiskGeomDet*>> HGCDiskGeomDetVector;

#endif