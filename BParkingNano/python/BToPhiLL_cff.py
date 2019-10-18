import FWCore.ParameterSet.Config as cms
from PhysicsTools.BParkingNano.common_cff import *


########## inputs preparation ################
electronPairsForPhiEE = cms.EDProducer(
    'DiElectronBuilder',
    src = cms.InputTag('electronsForAnalysis', 'SelectedElectrons'),
    transientTracksSrc = cms.InputTag('electronsForAnalysis', 'SelectedTransientElectrons'),
    lep1Selection = cms.string('pt > 1.5 && userFloat("unBiased") >= 3'),
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
       'PhiBuilder',
        pfcands= cms.InputTag('tracksBPark', 'SelectedTracks'),
        transientTracks= cms.InputTag('tracksBPark', 'SelectedTransientTracks'),
        trk1Selection = cms.string('pt > 1.5 && abs(eta)<2.4'), #need optimization   
        trk2Selection = cms.string('pt > 1.0 && abs(eta)<2.4'), #need optimization
        preVtxSelection = cms.string('abs(userCand("trk1").vz - userCand("trk2").vz)<1.0' 
        ' &&  pt()>2.0 && (mass() < 1.08 && mass() > 0.98)'
        ),
        postVtxSelection = cms.string('userFloat("sv_prob") > 1.e-5'
        ' && (userFloat("fitted_mass")<1.08 && userFloat("fitted_mass")>0.98)'
)
)



########################### B-> phi ll ##########################
BToPhiMuMu = cms.EDProducer(
    'BToPhiLLBuilder',
    dileptons = cms.InputTag('muonPairsForPhiMuMu'),
    leptonTransientTracks = muonPairsForPhiMuMu.transientTracksSrc,
    phis = cms.InputTag('PhiToKK'),
    phisTransientTracks = cms.InputTag('tracksBPark', 'SelectedTransientTracks'),
    
    beamSpot = cms.InputTag("offlineBeamSpot"),
    preVtxSelection = cms.string(
        'pt > 3. && userFloat("min_dr") > 0.03'
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
    
    beamSpot = cms.InputTag("offlineBeamSpot"),
    preVtxSelection = cms.string(
        'pt > 3. && userFloat("min_dr") > 0.03'
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
        # Mll
        mll_raw = Var('userCand("dilepton").mass()', float),
        mll_llfit = Var('userCand("dilepton").userFloat("fitted_mass")', float),
        mll_fullfit = ufloat('mll_fullfit'),     
        # phi fitted in b0 vertex
        mphi_fullfit = ufloat('mphi_fullfit'),
        ptphi_fullfit = ufloat('ptphi_fullfit'),
        etaphi_fullfit = ufloat('etaphi_fullfit'),
        phiphi_fullfit = ufloat('phiphi_fullfit'),
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

