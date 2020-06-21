import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *



trgTable = cms.EDProducer( "TrgBitTableProducer",
                          hltresults = cms.InputTag("TriggerResults::HLT"),
                          l1results  = cms.InputTag("gtStage2Digis::RECO"),
                          #add interesting paths
                          paths      = cms.vstring(
                                             "HLT_Mu7_IP4",
                                             "HLT_Mu8_IP6",
                                             "HLT_Mu8_IP5",
                                             "HLT_Mu8_IP3",
                                             "HLT_Mu8p5_IP3p5",
                                             "HLT_Mu9_IP6",
                                             "HLT_Mu9_IP5",
                                             "HLT_Mu9_IP4",    
                                             "HLT_Mu10p5_IP3p5",
                                             "HLT_Mu12_IP6",
                                             "HLT_DoubleMu4_3_Bs",
                                             "HLT_DoubleMu4_3_Jpsi_Displaced",
                                             "HLT_DoubleMu4_JpsiTrk_Displaced",
                                             "HLT_DoubleMu4_LowMassNonResonantTrk_Displaced",
                                             "HLT_DoubleMu4_PsiPrimeTrk_Displaced",
                                             "HLT_DoubleMu4_Jpsi_Displaced",
                                             "HLT_DoubleMu4_JpsiTrkTrk_Displaced",
                                             "HLT_Dimuon0_Jpsi_Muon",
                                             "HLT_Dimuon0er16_Jpsi_NoOS_NoVertexing",
                                             "HLT_Dimuon0er16_Jpsi_NoVertexing",
                                             "HLT_Dimuon10_Jpsi_Barrel",
                                             "HLT_Dimuon13_PsiPrime",
                                             "HLT_Dimuon16_Jpsi",
                                             "HLT_Dimuon20_Jpsi",
                                             "HLT_Dimuon6_Jpsi_NoVertexing",
                                             "HLT_Dimuon8_PsiPrime_Barrel",
                                             "HLT_DoubleMu3_Trk_Tau3mu",
                                             "HLT_Mu7p5_L2Mu2_Jpsi",
                                             "HLT_Mu7p5_Track3p5_Jpsi",
                                             "HLT_Mu7p5_Track7_Jpsi",
                                             "HLT_QuadMuon0_Dimuon0_Jpsi"
                                              ),
                           #add interesting seeds
                           seeds     = cms.vstring(
                                             "L1_SingleMu7er1p5",
                                             "L1_SingleMu8er1p5",
                                             "L1_SingleMu9er1p5",
                                             "L1_SingleMu10er1p5",
                                             "L1_SingleMu12er1p5",
                                             "L1_SingleMu22"
                                              ),
                            
)

trgTables = cms.Sequence(trgTable)



