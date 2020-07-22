import FWCore.ParameterSet.Config as cms
from PhysicsTools.BParkingNano.common_cff import *


########## inputs preparation ################
electronPairsForPhiEE = cms.EDProducer(
    'DiElectronBuilder',
    src = cms.InputTag('electronsForAnalysis', 'SelectedElectrons'),
    transientTracksSrc = cms.InputTag('electronsForAnalysis', 'SelectedTransientElectrons'),
    lep1Selection = cms.string('pt > 1.5 && userFloat("mvaId") >= 3'),
    lep2Selection = cms.string(''),
    preVtxSelection = cms.string(
        'abs(userCand("l1").vz - userCand("l2").vz) <= 1. && mass() < 5 '
        '&& mass() > 0 && charge() == 0 && userFloat("lep_deltaR") > 0.03'
    ),
    postVtxSelection = cms.string('userFloat("sv_chi2") < 998 && userFloat("sv_prob") > 1.e-5'),
)

muonPairsForPhiMuMu = cms.EDProducer(
    'DiMuonBuilder',
    src = cms.InputTag('muonTrgSelector', 'SelectedMuons'),
    transientTracksSrc = cms.InputTag('muonTrgSelector', 'SelectedTransientMuons'),
    lep1Selection = cms.string('pt > 1.5'),
    lep2Selection = cms.string(''),
    preVtxSelection = cms.string('abs(userCand("l1").vz - userCand("l2").vz) <= 1. && mass() < 5 '
                                 '&& mass() > 0 && charge() == 0 && userFloat("lep_deltaR") > 0.03'),
    postVtxSelection = electronPairsForPhiEE.postVtxSelection,
)

PhiToKK = cms.EDProducer(
       'DiTrackBuilder',
        pfcands= cms.InputTag('tracksBPark', 'SelectedTracks'),
        transientTracks= cms.InputTag('tracksBPark', 'SelectedTransientTracks'),
        trk1Mass = cms.double(K_MASS),
        trk2Mass = cms.double(K_MASS),
        trk1Selection = cms.string('pt > 0.7 && abs(eta)<2.5'), #need optimization   
        trk2Selection = cms.string('pt > 0.7 && abs(eta)<2.5'), #need optimization
        preVtxSelection = cms.string('abs(userCand("trk1").vz - userCand("trk2").vz)<1.0' 
        ' && pt()>0.0 && (mass() < 1.08 && mass() > 0.96)'
        ),
        postVtxSelection = cms.string('userFloat("sv_chi2") < 998 && userFloat("sv_prob") > 1.e-5'
        ' && (userFloat("fitted_mass")<1.06 && userFloat("fitted_mass")>0.98)'
)
)



########################### B-> phi ll ##########################
BToPhiMuMu = cms.EDProducer(
    'BToPhiLLBuilder',
    dileptons = cms.InputTag('muonPairsForPhiMuMu'),
    leptonTransientTracks = muonPairsForPhiMuMu.transientTracksSrc,
    phis = cms.InputTag('PhiToKK'),
    phisTransientTracks = cms.InputTag('tracksBPark', 'SelectedTransientTracks'),
    tracks = cms.InputTag("packedPFCandidates"),
    lostTracks = cms.InputTag("lostTracks"),
    isoTracksSelection = cms.string('pt > 0.7 && abs(eta)<2.5'),
    
    beamSpot = cms.InputTag("offlineBeamSpot"),
    preVtxSelection = cms.string(
        'pt > 3. && userFloat("min_dr") > 0.0'
        '&& (mass < 7. && mass > 4.) '
        ),
    postVtxSelection = cms.string(
        'userFloat("sv_prob") > 0.001 '
        '&& userFloat("fitted_cos_theta_2D") >= 0'
        '&& (userFloat("fitted_mass") > 4.5 && userFloat("fitted_mass") < 6.)'
    )
)

BToPhiEE = cms.EDProducer(
    'BToPhiLLBuilder',
    dileptons = cms.InputTag('electronPairsForPhiEE'),
    leptonTransientTracks = electronPairsForPhiEE.transientTracksSrc,
    phis = cms.InputTag('PhiToKK'),
    phisTransientTracks = cms.InputTag('tracksBPark', 'SelectedTransientTracks'),
    tracks = cms.InputTag("packedPFCandidates"),
    lostTracks = cms.InputTag("lostTracks"),
    isoTracksSelection = BToPhiMuMu.isoTracksSelection,
    
    beamSpot = cms.InputTag("offlineBeamSpot"),
    preVtxSelection = cms.string(
        'pt > 3. && userFloat("min_dr") > 0.0'
        '&& (mass < 7. && mass > 4.) '
        ),
    postVtxSelection = cms.string(
        'userFloat("sv_prob") > 0.001 '
        '&& userFloat("fitted_cos_theta_2D") >= 0'
        '&& (userFloat("fitted_mass") > 4.5 && userFloat("fitted_mass") < 6.)'
    )
)


################################### Tables #####################################

PhiToKKTable = cms.EDProducer(
   'SimpleCompositeCandidateFlatTableProducer',
    src = cms.InputTag("PhiToKK"),
    cut = cms.string(""),
    name = cms.string("Phi"),
    doc = cms.string("Phi Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=cms.PSet(
      CandVars,
      fitted_mass = ufloat('fitted_mass'),
      fitted_pt = ufloat('fitted_pt'),
      fitted_eta = ufloat('fitted_eta'),
      fitted_phi = ufloat('fitted_phi'),
      svprob = ufloat('sv_prob'),         
      trk_deltaR = ufloat('trk_deltaR'),
      trk1_idx = uint('trk1_idx'),
      trk2_idx = uint('trk2_idx')
    )
)


BToPhiEETable = cms.EDProducer(
    'SimpleCompositeCandidateFlatTableProducer',
    src = cms.InputTag("BToPhiEE"),
    cut = cms.string(""),
    name = cms.string("BToPhiEE"),
    doc = cms.string("BToPhiEE Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=cms.PSet(
        # pre-fit quantities
        CandVars,
        l1_idx = uint('l1_idx'),
        l2_idx = uint('l2_idx'),
        trk1_idx = uint('trk1_idx'),
        trk2_idx = uint('trk2_idx'),
        phi_idx = uint('phi_idx'),
        min_dr = ufloat('min_dr'),
        max_dr = ufloat('max_dr'),
        # fit and vtx info
        chi2 = ufloat('sv_chi2'),
        svprob = ufloat('sv_prob'),
        l_xy = ufloat('l_xy'),
        l_xy_unc = ufloat('l_xy_unc'),
        vtx_x = ufloat('vtx_x'),
        vtx_y = ufloat('vtx_y'),
        vtx_z = ufloat('vtx_z'),
        vtx_ex = ufloat('vtx_ex'), ## only saving diagonal elements of the cov matrix
        vtx_ey = ufloat('vtx_ey'),
        vtx_ez = ufloat('vtx_ez'),
        # Mll
        mll_raw = Var('userCand("dilepton").mass()', float),
        mll_llfit = Var('userCand("dilepton").userFloat("fitted_mass")', float),
        mll_fullfit = ufloat('fitted_mll'),     
        # phi fitted in b0 vertex
        fit_phi_mass = ufloat('fitted_phi_mass'),
        fit_phi_pt = ufloat('fitted_phi_pt'),
        fit_phi_eta = ufloat('fitted_phi_eta'),
        fit_phi_phi = ufloat('fitted_phi_phi'),
        # Cos(theta)
        cos2D = ufloat('cos_theta_2D'),
        fit_cos2D = ufloat('fitted_cos_theta_2D'),
        # post-fit momentum
        fit_mass = ufloat('fitted_mass'),
        fit_massErr = ufloat('fitted_massErr'),
        fit_pt = ufloat('fitted_pt'),
        fit_eta = ufloat('fitted_eta'),
        fit_phi = ufloat('fitted_phi'),
        # post-fit tracks/leptons
        #l1
        fit_l1_pt  = ufloat('fitted_l1_pt'),
        fit_l1_eta = ufloat('fitted_l1_eta'),
        fit_l1_phi = ufloat('fitted_l1_phi'),
        #l2
        fit_l2_pt  = ufloat('fitted_l2_pt'),
        fit_l2_eta = ufloat('fitted_l2_eta'),
        fit_l2_phi = ufloat('fitted_l2_phi'),
        #trk1
        fit_trk1_pt  = ufloat('fitted_trk1_pt'),
        fit_trk1_eta = ufloat('fitted_trk1_eta'),
        fit_trk1_phi = ufloat('fitted_trk1_phi'),
        #trk2
        fit_trk2_pt  = ufloat('fitted_trk2_pt'),
        fit_trk2_eta = ufloat('fitted_trk2_eta'),
        fit_trk2_phi = ufloat('fitted_trk2_phi'),
        # isolation 
        l1_iso03 = ufloat('l1_iso03'),
        l1_iso04 = ufloat('l1_iso04'),
        l2_iso03 = ufloat('l2_iso03'),
        l2_iso04 = ufloat('l2_iso04'),
        trk1_iso03 = ufloat('trk1_iso03'),
        trk1_iso04 = ufloat('trk1_iso04'),
        trk2_iso03 = ufloat('trk2_iso03'),
        trk2_iso04 = ufloat('trk2_iso04'),
        b_iso03  = ufloat('b_iso03'),
        b_iso04  = ufloat('b_iso04'),
    )
)

BToPhiMuMuTable = BToPhiEETable.clone(
    src = cms.InputTag("BToPhiMuMu"),
    name = cms.string("BToPhiMuMu"),
    doc = cms.string("BToPhiMuMu Variables")
)

CountBToPhiEE = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("BToPhiEE")
)    
CountBToPhiMuMu = CountBToPhiEE.clone(
    minNumber = cms.uint32(1),
    src = cms.InputTag("BToPhiMuMu")
)


########################### Sequencies  ############################

PhiToKKSequence = cms.Sequence(  PhiToKK  )

BToPhiMuMuSequence = cms.Sequence(
    (muonPairsForPhiMuMu *BToPhiMuMu )
)


BToPhiEESequence = cms.Sequence(
    (electronPairsForPhiEE *BToPhiEE )
)


BToPhiLLSequence = cms.Sequence(
    ( (muonPairsForPhiMuMu *BToPhiMuMu)
     +(electronPairsForPhiEE *BToPhiEE) )   
)


BToPhiLLTables = cms.Sequence( BToPhiEETable + BToPhiMuMuTable )

