import FWCore.ParameterSet.Config as cms
from RecoEgamma.ElectronIdentification.Identification.mvaElectronID_tools import *
from os import path

# Egamma presentation on this ID:
# https://indico.cern.ch/event/732971/contributions/3022864/attachments/1658765/2656595/180530_egamma.pdf

mvaTag = "BParkRetrain"

weightFileDir = "PhysicsTools/BParkingNano/data/PFRetrainWeightFiles"

mvaWeightFiles = [
     path.join(weightFileDir, "BParkRetrain_LowPt_unbiased.xml.gz"),
     path.join(weightFileDir, "BParkRetrain_HighPt_unbiased.xml.gz"),
     ]

EleMVA_2CategoriesCuts = [
    "pt < 5.",
    "pt >= 5.",
    ]

mvaEleID_BParkRetrain_producer_config = cms.PSet(
    mvaName             = cms.string(mvaClassName),
    mvaTag              = cms.string(mvaTag),
    nCategories         = cms.int32(2),
    categoryCuts        = cms.vstring(*EleMVA_2CategoriesCuts),
    weightFileNames     = cms.vstring(*mvaWeightFiles),
    variableDefinition  = cms.string(mvaVariablesFile)
    )

