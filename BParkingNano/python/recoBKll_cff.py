import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *

process.reconstructB = cms.EDFilter(
        "BKllProducer",
         muons = cms.InputTag("muonsBPark:SelectedMuons"),
         electrons = cms.InputTag("electronsBPark:SelectedElectrons"),
         tracks = cms.InputTag("tracksBPark:SelectedTracks"),
         beamSpot= cms.InputTag("offlineBeamSpot"),
         RunParameters = cms.PSet(
           SkipNoRecoBEvt = cms.bool(False),
           pdgId = cms.int32(13),
           name = cms.string("BKmumu"),
           ptLep1Cut = cms.double(0),
           mllMax= cms.double(5.0),
           mllMin= cms.double(0),
           MBMin = cms.double(4.5),
           MBMax = cms.double(6),
           drTrkLepCut = cms.double(0.01)
         )
)




bToKMuMuTable = cms.EDProducer(
    'SimpleCompositeCandidateFlatTableProducer',
    src = cms.InputTag("reconstructB:BKmumu"),
    cut = cms.string(""),
    name = cms.string("BToKMuMu"),
    doc = cms.string("BToKMuMu Variable"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=cms.PSet(
        CandVars,
    )
)


recoBseq= cms.Sequence (
   reconstructB
  +bToKMuMuTable
)




