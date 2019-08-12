import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *



trgTable = cms.EDProducer( "TrgBitTableProducer",
                          triggerresults= cms.InputTag("TriggerResults::HLT"),
                           paths= cms.vstring("HLT_Mu9_IP6",    "HLT_Mu9_IP5", 
                                              "HLT_Mu7_IP4",    "HLT_Mu8_IP3",
                                              "HLT_Mu12_IP6"
                                              ),
)

trgTable = cms.Sequence(trgTable)



