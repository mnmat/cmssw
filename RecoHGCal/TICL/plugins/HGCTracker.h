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

		const HGCDiskLayer* nextDisk(const HGCDiskLayer * from, 
                                        PropagationDirection direction, 
      									bool isSilicon) const;

		const HGCDiskLayer* disk(int zside, 
									int disk) const { 
					return (zside > 0 ? diskLayerPos_ : diskLayerNeg_).at(disk); };

        const HGCDiskLayer* firstDisk(int zside) const { return (zside > 0 ? diskLayerPos_ : diskLayerNeg_).front(); }
        const HGCDiskLayer* lastDisk(int zside) const { return (zside > 0 ? diskLayerPos_ : diskLayerNeg_).back(); }
        const HGCDiskLayer* firstDisk(int zside, PropagationDirection direction) const { return direction == alongMomentum ? firstDisk(zside) : lastDisk(zside); }
        const HGCDiskLayer* lastDisk(int zside, PropagationDirection direction) const { return direction == alongMomentum ? lastDisk(zside) : firstDisk(zside); }
	
	private:
		// Member Functions


		void makeDisks(int subdet, const CaloGeometry* geom, std::vector<HGCDiskGeomDet*>& disksPos, std::vector<HGCDiskGeomDet*>& diskssNeg);
		void makeDiskLayers(std::vector<HGCDiskGeomDet*>disksSc, std::vector<HGCDiskGeomDet*> disksSi);

		void addDiskLayer(HGCDiskGeomDet *disk) { 
      		(disk->zside() > 0 ? diskLayerPos_ : diskLayerNeg_).push_back(new HGCDiskLayer(nullptr, disk));
      	}

		void addDiskLayer(HGCDiskGeomDet *silicon, HGCDiskGeomDet *scintillator) { 
      		(silicon->zside() > 0 ? diskLayerPos_ : diskLayerNeg_).push_back(new HGCDiskLayer(scintillator, silicon));
      	}
		// Member Variables

		std::vector<HGCDiskGeomDet *> disksPos_, disksNeg_;
		std::vector<HGCDiskLayer*> diskLayerPos_, diskLayerNeg_;

		hgcal::RecHitTools rhtools_;
		const CaloGeometry* geom_;

		const std::vector<double> radlen_, xi_;
		int lastLayer = 47; // use rhtools to get values
		int offset = 34; // Use rhtools to get values
};

#endif
