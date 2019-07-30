from PhysicsTools.NanoAOD.common_cff import *
import FWCore.ParameterSet.Config as cms

RecoTrgMuonCand=cms.EDProducer("TriggeringMuonProducer",
                               muonCollection = cms.InputTag("slimmedMuons"), #same collection as in NanoAOD
                               bits = cms.InputTag("TriggerResults","","HLT"),
                               prescales = cms.InputTag("patTrigger"),
                               objects = cms.InputTag("slimmedPatTrigger"),
                               vertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
                               maxdR_matching = cms.double(0.01)
                               )


#RecoTrgMuonCandTable=cms.EDProducer("SimpleMuonFlatTableProducer",
RecoTrgMuonCandTable=cms.EDProducer("SimpleCandidateFlatTableProducer",
                                    src=cms.InputTag("RecoTrgMuonCand"),
                                    cut=cms.string(""),
                                    name=cms.string("RecoTrgMuonCand"),
                                    doc=cms.string("reco muon matched to triggering muon"),
                                    singleton=cms.bool(False),
                                    extension=cms.bool(False),
                                    variables=cms.PSet(recoMuonIndex = Var("userInt('recoMuonIndex')", int,doc="index in original reco muon collection"),
                                                       trgMuonIndex = Var("userInt('trgMuonIndex')", int,doc="index in trigger muon collection")
                                                       #these following are already in the full muon collection
                                                       #pt = Var("daughter(0).pt()",float,doc="pt"),
                                                       #eta = Var("daughter(0).eta()",float,doc="eta"),
                                                       #phi = Var("daughter(0).phi()",float,doc="phi"),
                                                       #vz = Var("daughter(0).vz()", float, doc="z coordinate of vertex position, in cm"),
                                                   ),
                                    )


RecoTrgMuonCandSequence=cms.Sequence(RecoTrgMuonCand)
RecoTrgMuonCandTables=cms.Sequence(RecoTrgMuonCandTable)
