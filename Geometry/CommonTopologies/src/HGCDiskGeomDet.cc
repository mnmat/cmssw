#include "Geometry/CommonTopologies/interface/HGCDiskGeomDet.h"
#include "DataFormats/GeometrySurface/interface/MediumProperties.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"

HGCDiskGeomDet::HGCDiskGeomDet(int subdet, int zside, int layer, float z, float rmin, float rmax, float radlen, float xi) :
            GeomDet( Disk::build(Disk::PositionType(0,0,z), Disk::RotationType(), SimpleDiskBounds(rmin, rmax, -20, 20)).get() ),
            subdet_(subdet), zside_(zside), layer_(layer), rmin_(rmin), rmax_(rmax) 
{
    if (radlen > 0) {
        (const_cast<Plane &>(surface())).setMediumProperties(MediumProperties(radlen,xi));
    }
}

#include "FWCore/Utilities/interface/typelookup.h"
TYPELOOKUP_DATA_REG(HGCDiskGeomDet);
TYPELOOKUP_DATA_REG(HGCDiskGeomDetVector);