import FWCore.ParameterSet.Config as cms
from PhysicsTools.BParkingNano.common_cff import *

electronPairsForKee = cms.EDProducer(
    'DiElectronBuilder',
    src = cms.InputTag('electronsForAnalysis', 'SelectedElectrons'),
    transientTracksSrc = cms.InputTag('electronsForAnalysis', 'SelectedTransientElectrons'),
    lep1Selection = cms.string('pt > 1.5 && userFloat("unBiased") >= 3'),
    lep2Selection = cms.string(''),
    min_dr = cms.double(0.03), ##for the moment than we will see
    preVtxSelection = cms.string(
        'abs(userCand("l1").vz - userCand("l2").vz) <= 1. && mass() < 5 '
        '&& mass() > 0 && charge() == 0 && userFloat("lep_deltaR") > 0.03'
    ),
    postVtxSelection = cms.string('userFloat("sv_chi2") < 998 && userFloat("sv_prob") > 1.e-5'),
    isMuonTrack = cms.bool(False)
)

BToKee = cms.EDProducer(
    'BToKLLBuilder',
    dileptons = cms.InputTag('electronPairsForKee'),
    leptonTransientTracks = electronPairsForKee.transientTracksSrc,
    kaons = cms.InputTag('tracksBPark', 'SelectedTracks'),
    kaonsTransientTracks = cms.InputTag('tracksBPark', 'SelectedTransientTracks'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    kaonSelection = cms.string(''),
    preVtxSelection = cms.string(
        'pt > 3. && userFloat("min_dr") > 0.03 '
        '&& mass < 7. && mass > 4.'
        ),
    postVtxSelection = cms.string(
        'userInt("sv_OK") == 1 && userFloat("sv_prob") > 0.001 '
        '&& userFloat("fitted_cos_theta_2D") >= 0 '
        '&& userFloat("fitted_mass") > 4.5 && userFloat("fitted_mass") < 6.'
    )
)

muonPairsForKmumu = cms.EDProducer(
    'DiMuonBuilder',
    src = cms.InputTag('muonTrgSelector', 'SelectedMuons'),
    transientTracksSrc = cms.InputTag('muonTrgSelector', 'SelectedTransientMuons'),
    lep1Selection = cms.string('pt > 1.5'),
    lep2Selection = cms.string(''),
    min_dr = cms.double(0.01),
    preVtxSelection = cms.string('abs(userCand("l1").vz - userCand("l2").vz) <= 1. && mass() < 5 '
                                 '&& mass() > 0 && charge() == 0 && userFloat("lep_deltaR") > 0.03'),
    postVtxSelection = electronPairsForKee.postVtxSelection,
    isMuonTrack = cms.bool(False)
)

muonTrackPairsForKmumu = cms.EDProducer(
    'DiMuonTrackBuilder',
    src = cms.InputTag("muonsTracksBPark:SelectedMuonsTracks"),
    transientTracksSrc = cms.InputTag('muonsTracksBPark', 'SelectedTransientMuonsTracks'),
    lep1Selection = cms.string('pt > 1.5'),
    lep2Selection = cms.string(''),
    min_dr = cms.double(0.01),
    preVtxSelection = cms.string('abs(userCand("l1").vz - userCand("l2").vz) <= 1. && mass() < 5 '
                                 '&& mass() > 0 && charge() == 0'),
    postVtxSelection = electronPairsForKee.postVtxSelection,
    isMuonTrack = cms.bool(True)
)


BToKmumu = cms.EDProducer(
    'BToKLLBuilder',
    #dileptons = cms.InputTag('muonPairsForKmumu'),
    #leptonTransientTracks = muonPairsForKmumu.transientTracksSrc,
    dileptons = cms.InputTag('muonTrackPairsForKmumu'),
    leptonTransientTracks = muonTrackPairsForKmumu.transientTracksSrc,
    kaons = BToKee.kaons,
    kaonsTransientTracks = BToKee.kaonsTransientTracks,
    beamSpot = cms.InputTag("offlineBeamSpot"),
    kaonSelection = cms.string(''),
    # This in principle can be different between electrons and muons
    preVtxSelection = cms.string(
        'pt > 3. && userFloat("min_dr") > 0.03'
        '&& mass < 7. && mass > 4.'
        ),
    postVtxSelection = cms.string(
        'userInt("sv_OK") == 1 && userFloat("sv_prob") > 0.001 '
        '&& userFloat("fitted_cos_theta_2D") >= 0'
        '&& userFloat("fitted_mass") > 4.5 && userFloat("fitted_mass") < 6.'
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
        # pre-fit quantities
        CandVars,
        l1Idx = uint('l1_idx'),
        l2Idx = uint('l2_idx'),
        l1OrigIdx = uint('l1_original_idx'),
        l2OrigIdx = uint('l2_original_idx'),
        kIdx = uint('k_idx'),
        l1_isPF = uint('l1_isPF'),
        l1_isTrack = uint('l1_isTrack'),
        l2_isPF = uint('l2_isPF'),
        l2_isTrack = uint('l2_isTrack'),
        minDR = ufloat('min_dr'),
        maxDR = ufloat('max_dr'),
        # fit and vtx info
        chi2 = ufloat('sv_chi2'),
        svprob = ufloat('sv_prob'),
        l_xy = ufloat('l_xy'),
        l_xy_unc = ufloat('l_xy_unc'),
        # Mll
        mll_raw = Var('userCand("dilepton").mass()', float),
        mll_llfit = Var('userCand("dilepton").userFloat("fitted_mass")', float), # this might not work
        mll_fullfit = ufloat('fitted_mll'),
        # Cos(theta)
        cos2D = ufloat('cos_theta_2D'),
        fit_cos2D = ufloat('fitted_cos_theta_2D'),
        # post-fit momentum
        fit_mass = ufloat('fitted_mass'),
        fit_pt = ufloat('fitted_pt'),
        fit_eta = ufloat('fitted_eta'),
        fit_phi = ufloat('fitted_phi'),
    )
)

BToKmumuTable = BToKeeTable.clone(
    src = cms.InputTag("BToKmumu"),
    name = cms.string("BToKMuMu"),
    doc = cms.string("BToKMuMu Variable")
)


CountBToKee = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("BToKee")
)    
CountBToKmumu = CountBToKee.clone(
    minNumber = cms.uint32(1),
    src = cms.InputTag("BToKmumu")
)


BToKMuMuSequence = cms.Sequence(
   # (muonPairsForKmumu * BToKmumu)
    (muonTrackPairsForKmumu * BToKmumu)
)
BToKEESequence = cms.Sequence(
    (electronPairsForKee * BToKee)
)

BToKLLSequence = cms.Sequence(
    (electronPairsForKee * BToKee) +
    (muonPairsForKmumu * BToKmumu)
)
BToKLLTables = cms.Sequence(BToKeeTable + BToKmumuTable)

