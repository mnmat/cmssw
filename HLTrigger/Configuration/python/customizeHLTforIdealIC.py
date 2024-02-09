import FWCore.ParameterSet.Config as cms

def customizeHLTforIdealIC(process):
  """Customisation to apply the ideal intercalibration constants to ECAL"""

  if hasattr(process, "GlobalTag"): 
    if not hasattr(process.GlobalTag, "toGet"):
      process.GlobalTag.toGet = cms.VPSet()
    process.GlobalTag.toGet += [
      cms.PSet(
        record = cms.string("EcalIntercalibConstantsRcd"),
        tag = cms.string("EcalIntercalibConstants_MC_Digi_2018"),
      )
    ]
  else:
    print(" customizeHLTforIdealIC -- the process.GlobalTag ESSource does not exist: no customisation applied.")
  return process 
