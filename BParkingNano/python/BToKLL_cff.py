import FWCore.ParameterSet.Config as cms
from PhysicsTools.BParkingNano.common_cff import *

electronPairsForKee = cms.EDProducer(
    'DiElectronBuilder',
    src = cms.InputTag('electronsForAnalysis', 'SelectedElectrons'),
    transientTracksSrc = cms.InputTag('electronsForAnalysis', 'FIXME'),
    lep1Selection = cms.string('pt > 1.5 && (userInt("isPF") == 1 || electronID("unbiased_seed") >= 3 )'),
    lep2Selection = cms.string('pt > 0.5 && (userInt("isPF") == 1 || electronID("unbiased_seed") >= -4)'),
    preVtxSelection = cms.string(
        'abs(userCand("l1").vz - userCand("l2").vz) <= 1. && mass() < 5 '
        '&& mass() > 0 && charge() == 0 && userFloat("lep_deltaR") > 0.03'
    ),
    postVtxSelection = cms.string('userFloat("sv_chi2") < 998 && userFloat("sv_prob") > 0'),
)

BToKee = cms.EDProducer(
    'BToKLLBuilder',
    dileptons = cms.InputTag('electronPairsForKee'),
    leptonTransientTracks = electronPairsForKee.transientTracksSrc,
    kaons = cms.InputTag('tracksBPark', 'SelectedTracks'),
    kaonsTransientTracks = cms.InputTag('tracksBPark', 'FIXME'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    kaonSelection = cms.string(''),
    preVtxSelection = cms.string(
        'pt > 3. && userFloat("min_dr") > 0.01 '
        '&& mass < 6 && mass > 4.5'
        ),
    postVtxSelection = cms.string(
        'userInt("sv_OK") == 1 && userFloat("sv_prob") > 0.00000001 '
        '&& userFloat("cos_theta_2D") >= 0'
    )
)

muonPairsForKmumu = cms.EDProducer(
    'DiMuonBuilder',
    src = cms.InputTag('muonTrgSelector', 'SelectedMuons'),
    transientTracksSrc = cms.InputTag('muonTrgSelector', 'FIXME'),
    lep1Selection = cms.string(''), #FIXME
    lep2Selection = cms.string(''),
    preVtxSelection = electronPairsForKee.preVtxSelection,
    postVtxSelection = electronPairsForKee.postVtxSelection,
)

BToKmumu = cms.EDProducer(
    'BToKLLBuilder',
    dileptons = cms.InputTag('muonPairsForKmumu'),
    leptonTransientTracks = muonPairsForKmumu.transientTracksSrc,
    kaons = cms.InputTag('tracksBPark', 'SelectedTracks'),
    kaonsTransientTracks = cms.InputTag('tracksBPark', 'FIXME'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    kaonSelection = cms.string(''),
    # This in principle can be different between electrons and muons
    preVtxSelection = cms.string(
        'pt > 3. && userFloat("min_dr") > 0.01 '
        '&& mass < 6 && mass > 4.5'
        ),
    postVtxSelection = cms.string(
        'userInt("sv_OK") == 1 && userFloat("sv_prob") > 0.00000001 '
        '&& userFloat("cos_theta_2D") >= 0'
    )
)

BToKeeTable = cms.EDProducer(
    'SimpleCompositeCandidateFlatTableProducer',
    src = cms.InputTag("BToKee"),
    cut = cms.string(""),
    name = cms.string("BToKEE"),
    doc = cms.string("BToKEE Variable"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=cms.PSet(
        CandVars,
        el1Idx = uint('l1_idx'),
        el2Idx = uint('l2_idx'),
        kIdx = uint('k_idx'),
        chi2 = ufloat('sv_chi2'),
        svprob = ufloat('sv_prob'),
        mll_raw = Var('userCand("dilepton").mass()', float),
        mll_llfit = Var('userCand("dilepton").userFloat("fitted_mass")', float), # this might not work
        mll_fullfit = ufloat('fitted_mll'),
        cos2D = ufloat('cos_theta_2D'),
        fit_mass = ufloat('fitted_mass'),
        l_xy = ufloat('l_xy'),
        l_xy_unc = ufloat('l_xy_unc'),
    )
)

BToKmumuTable = BToKeeTable.clone(
    src = cms.InputTag("BToKmumu"),
    name = cms.string("BToKMuMu"),
    doc = cms.string("BToKMuMu Variable")
)

BToKLLSequence = cms.Sequence(
    (electronPairsForKee * BToKee) +
    (muonPairsForKmumu * BToKmumu)
)
BToKLLTables = cms.Sequence(BToKeeTable + BToKmumuTable)

