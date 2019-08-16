import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *


reconstructB = cms.EDProducer(
        "BKllProducer",
         muons = cms.InputTag("muonTrgSelector:SelectedMuons"),
         electrons = cms.InputTag("electronsForAnalysis:SelectedElectrons"),
         tracks = cms.InputTag("tracksBPark:SelectedTracks"),
         beamSpot= cms.InputTag("offlineBeamSpot"),
         RunParameters = cms.PSet(
             SkipNoRecoBEvt = cms.bool(False),
           FinalLeptonId = cms.int32(13),
           label = cms.string("BKmumu"),
           ptLep1Cut = cms.double(1.5),
           bdtEl1Cut = cms.double(3.0),
           dzLepLepCut = cms.double(1.0),
           mllMax= cms.double(5.0),
           mllMin= cms.double(0),
           MBMin = cms.double(4.5),
           MBMax = cms.double(6),
           drTrkLepCut = cms.double(0.01),
           BPtCut = cms.double(3),
           BVtxProbCut = cms.double(0.001),
           BCosACut = cms.double(0)
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
        mu1Idx = Var('userInt("l1_idx")', int),
        mu2Idx = Var('userInt("l2_idx")', int),
        kIdx   = Var('userInt("k_idx" )', int),
        vtx_prob = Var('userFloat("vtx_prob" )', float),
        Lxy = Var('userFloat("Lxy" )', float),
        eLxy = Var('userFloat("eLxy" )', float),
        cosA = Var('userFloat("cosA" )', float),
        mll = Var('userFloat("mll" )', float, precision=10),
        l1refit_pt = Var('userFloat("l1refit_pt" )', float),
        l1refit_eta = Var('userFloat("l1refit_eta" )', float, precision=10),
        l1refit_phi = Var('userFloat("l1refit_phi" )', float, precision=10),
        l2refit_pt = Var('userFloat("l2refit_pt" )', float),
        l2refit_eta = Var('userFloat("l2refit_eta" )', float, precision=10),
        l2refit_phi = Var('userFloat("l2refit_phi" )', float, precision=10),
        Krefit_pt = Var('userFloat("Hrefit_pt" )', float),
        Krefit_eta = Var('userFloat("Hrefit_eta" )', float, precision=10),
        Krefit_phi = Var('userFloat("Hrefit_phi" )', float, precision=10),
    )
)


recoBSequence= cms.Sequence (
   reconstructB
)
recoBTableSequence= cms.Sequence(
   bToKMuMuTable
)




