#ifndef __RecoHGCal_TICL_HGCTracker_H__
#define __RecoHGCal_TICL_HGCTracker_H__

#include "Geometry/CommonTopologies/interface/HGCDiskGeomDet.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"


class HGCTracker{
	public:
		HGCTracker(const CaloGeometry *geom,
						 const std::vector<double> radlen_,
						 const std::vector<double> xi_);

		~HGCTracker();

		const HGCDiskGeomDet* nextDisk(const HGCDiskGeomDet * from, 
                                        PropagationDirection direction, 
      									bool isSilicon) const;
    	const HGCDiskGeomDet* switchDisk(const HGCDiskGeomDet * from) const;

		const HGCDiskGeomDet* disk(int zside, 
									int disk) const { 
					return (zside > 0 ? disksPos_ : disksNeg_).at(disk); };

        const HGCDiskGeomDet* firstDisk(int zside) const { return (zside > 0 ? disksPos_ : disksNeg_).front(); }
        const HGCDiskGeomDet* lastDisk(int zside) const { return (zside > 0 ? disksPos_ : disksNeg_).back(); }
        const HGCDiskGeomDet* firstDisk(int zside, PropagationDirection direction) const { return direction == alongMomentum ? firstDisk(zside) : lastDisk(zside); }
        const HGCDiskGeomDet* lastDisk(int zside, PropagationDirection direction) const { return direction == alongMomentum ? lastDisk(zside) : firstDisk(zside); }
	
	private:
		// Member Functions

		void makeDisks(int subdet, const CaloGeometry* geom);
	 	void addDisk(HGCDiskGeomDet *disk) { 
      		(disk->zside() > 0 ? disksPos_ : disksNeg_).push_back(disk);
      	}
		// Member Variables

		std::vector<HGCDiskGeomDet *> disksPos_, disksNeg_;
		hgcal::RecHitTools rhtools_;
		const CaloGeometry* geom_;

		const std::vector<double> radlen_, xi_;
};

#endif
