import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.triggerObjects_cff import *


triggerObjectBParkTable = cms.EDProducer("TriggerObjectTableBParkProducer",
    name= cms.string("TrigObj"),
    src = cms.InputTag("unpackedPatTrigger"),
    l1Muon = cms.InputTag("gmtStage2Digis","Muon"),
    selections = cms.VPSet(
        cms.PSet(
            name = cms.string("Muon"),
            id = cms.int32(13),
            sel = cms.string("type(83) && pt > 5 && coll('hltIterL3MuonCandidates')"), 
            l1seed = cms.string("type(-81)"), l1deltaR = cms.double(0.5),
            l2seed = cms.string("type(83) && coll('hltL2MuonCandidates')"),  l2deltaR = cms.double(0.3),
            qualityBits = cms.string("filter('hltL3fL1s*Park*')"), qualityBitsDoc = cms.string("1 = Muon filters for BPH parking"),
        ),
    ),
)

triggerObjectBParkTables = cms.Sequence( unpackedPatTrigger + triggerObjectBParkTable )
