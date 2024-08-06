#ifndef Geometry_Record_CaloRecoGeometryRecord_h
#define Geometry_Record_CaloRecoGeometryRecord_h

/** \class CaloRecoGeometryRecord
 *
 *  Record to hold calo reconstruction geometries.
 *
 *  \author M. Matthewman - CERN
 */

#include "FWCore/Framework/interface/EventSetupRecordImplementation.h"
#include "FWCore/Framework/interface/DependentRecordImplementation.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "FWCore/Utilities/interface/mplVector.h"

class CaloRecoGeometryRecord
    : public edm::eventsetup::DependentRecordImplementation<CaloRecoGeometryRecord,
                                                            edm::mpl::Vector<CaloGeometryRecord> > {};

#endif
