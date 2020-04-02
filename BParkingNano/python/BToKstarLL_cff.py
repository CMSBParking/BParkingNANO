import FWCore.ParameterSet.Config as cms
from PhysicsTools.BParkingNano.common_cff import *


########## inputs preparation ################
electronPairsForKstarEE = cms.EDProducer(
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

muonPairsForKstarMuMu = cms.EDProducer(
    'DiMuonBuilder',
    src = cms.InputTag('muonTrgSelector', 'SelectedMuons'),
    transientTracksSrc = cms.InputTag('muonTrgSelector', 'SelectedTransientMuons'),
    lep1Selection = cms.string('pt > 1.5'),
    lep2Selection = cms.string(''),
    preVtxSelection = cms.string('abs(userCand("l1").vz - userCand("l2").vz) <= 1. && mass() < 5 '
                                 '&& mass() > 0 && charge() == 0 && userFloat("lep_deltaR") > 0.03'),
    postVtxSelection = electronPairsForKstarEE.postVtxSelection,
)

KstarToKPi = cms.EDProducer(
       'KstarBuilder',
        pfcands= cms.InputTag('tracksBPark', 'SelectedTracks'),
        transientTracks= cms.InputTag('tracksBPark', 'SelectedTransientTracks'),
        trk1Selection = cms.string('pt > 1.5 && abs(eta)<2.4'), #need optimization   
        trk2Selection = cms.string('pt > 1.0 && abs(eta)<2.4'), #need optimization
        preVtxSelection = cms.string('abs(userCand("trk1").vz - userCand("trk2").vz)<1.0' 
        ' &&  pt()>2.0 && ( (mass() < 1.042 && mass() > 0.742)'
        ' || (userFloat("barMass") < 1.042 && userFloat("barMass") > 0.742) ) '
        ),
        postVtxSelection = cms.string('userFloat("sv_prob") > 1.e-5'
        ' && (  (userFloat("fitted_mass")<1.042 && userFloat("fitted_mass")>0.742)'
        ' || (userFloat("fitted_barMass")<1.042 && userFloat("fitted_barMass")>0.742)  )'
)
)



########################### B-> K* ll ##########################
BToKstarMuMu = cms.EDProducer(
    'BToKstarLLBuilder',
    dileptons = cms.InputTag('muonPairsForKstarMuMu'),
    leptonTransientTracks = muonPairsForKstarMuMu.transientTracksSrc,
    kstars = cms.InputTag('KstarToKPi'),
    kstarsTransientTracks = cms.InputTag('tracksBPark', 'SelectedTransientTracks'),
    tracks = cms.InputTag("packedPFCandidates"),
    lostTracks = cms.InputTag("lostTracks"),
    isoTracksSelection = cms.string('pt > 0.7 && abs(eta)<2.5'),
    
    beamSpot = cms.InputTag("offlineBeamSpot"),
    preVtxSelection = cms.string(
        'pt > 3. && userFloat("min_dr") > 0.03'
        '&& ( (mass < 7. && mass > 4.) '
        '|| (userFloat("barMass")<7. && userFloat("barMass")>4.) )'
        ),
    postVtxSelection = cms.string(
        'userFloat("sv_prob") > 0.001 '
        '&& userFloat("fitted_cos_theta_2D") >= 0'
        '&& ( (userFloat("fitted_mass") > 4.5 && userFloat("fitted_mass") < 6.)'
        '|| (userFloat("fitted_barMass") > 4.5 && userFloat("fitted_barMass") < 6.)  )'
    )
)

BToKstarEE = cms.EDProducer(
    'BToKstarLLBuilder',
    dileptons = cms.InputTag('electronPairsForKstarEE'),
    leptonTransientTracks = electronPairsForKstarEE.transientTracksSrc,
    kstars = cms.InputTag('KstarToKPi'),
    kstarsTransientTracks = cms.InputTag('tracksBPark', 'SelectedTransientTracks'),
    tracks = cms.InputTag("packedPFCandidates"),
    lostTracks = cms.InputTag("lostTracks"),
    isoTracksSelection = BToKstarMuMu.isoTracksSelection,
    
    beamSpot = cms.InputTag("offlineBeamSpot"),
    preVtxSelection = cms.string(
        'pt > 3. && userFloat("min_dr") > 0.03'
        '&& ( (mass < 7. && mass > 4.) '
        '|| (userFloat("barMass")<7. && userFloat("barMass")>4.) )'
        ),
    postVtxSelection = cms.string(
        'userFloat("sv_prob") > 0.001 '
        '&& userFloat("fitted_cos_theta_2D") >= 0'
        '&& ( (userFloat("fitted_mass") > 4.5 && userFloat("fitted_mass") < 6.)'
        '|| (userFloat("fitted_barMass") > 4.5 && userFloat("fitted_barMass") < 6.)  )'
    )
)


################################### Tables #####################################

KstarToKPiTable = cms.EDProducer(
   'SimpleCompositeCandidateFlatTableProducer',
    src = cms.InputTag("KstarToKPi"),
    cut = cms.string(""),
    name = cms.string("Kstar"),
    doc = cms.string("Kstar Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=cms.PSet(
      CandVars,
      barMass = ufloat('barMass'),
      fitted_mass = ufloat('fitted_mass'),
      fitted_barMass = ufloat('fitted_barMass'),
      fitted_pt = ufloat('fitted_pt'),
      fitted_eta = ufloat('fitted_eta'),
      fitted_phi = ufloat('fitted_phi'),
      svprob = ufloat('sv_prob'),         
      trk_deltaR = ufloat('trk_deltaR'),
      trk1_idx = uint('trk1_idx'),
      trk2_idx = uint('trk2_idx')
    )
)


BToKstarEETable = cms.EDProducer(
    'SimpleCompositeCandidateFlatTableProducer',
    src = cms.InputTag("BToKstarEE"),
    cut = cms.string(""),
    name = cms.string("BToKsEE"),
    doc = cms.string("BToKstarEE Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=cms.PSet(
        # pre-fit quantities
        CandVars,
        l1_idx = uint('l1_idx'),
        l2_idx = uint('l2_idx'),
        trk1_idx = uint('trk1_idx'),
        trk2_idx = uint('trk2_idx'),
        kstar_idx = uint('kstar_idx'),
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
        mll_fullfit = ufloat('fitted_mll'),     
        # kstar fitted in b0 vertex
        fit_kstar_mass = ufloat('fitted_kstar_mass'),
        fit_kstar_pt = ufloat('fitted_kstar_pt'),
        fit_kstar_eta = ufloat('fitted_kstar_eta'),
        fit_kstar_phi = ufloat('fitted_kstar_phi'),
        # Cos(theta)
        cos2D = ufloat('cos_theta_2D'),
        fit_cos2D = ufloat('fitted_cos_theta_2D'),
        # post-fit momentum
        fit_mass = ufloat('fitted_mass'),
        fit_massErr = ufloat('fitted_massErr'),
        fit_pt = ufloat('fitted_pt'),
        fit_eta = ufloat('fitted_eta'),
        fit_phi = ufloat('fitted_phi'),
        # additional mass hypothesis
        barMass = ufloat ('barMass'),
        fit_barMass = ufloat('fitted_barMass'),
        fit_barKstar_mass = ufloat('barMasskstar_fullfit'),
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
        tk1_iso03 = ufloat('tk1_iso03'),
        tk1_iso04 = ufloat('tk1_iso04'),
        tk2_iso03 = ufloat('tk2_iso03'),
        tk2_iso04 = ufloat('tk2_iso04'),
        b_iso03  = ufloat('b_iso03'),
        b_iso04  = ufloat('b_iso04'),
    )
)

BToKstarMuMuTable = BToKstarEETable.clone(
    src = cms.InputTag("BToKstarMuMu"),
    name = cms.string("BToKsMuMu"),
    doc = cms.string("BToKstarMuMu Variables")
)

CountBToKstarEE = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("BToKstarEE")
)    
CountBToKstarMuMu = CountBToKstarEE.clone(
    minNumber = cms.uint32(1),
    src = cms.InputTag("BToKstarMuMu")
)


########################### Sequencies  ############################

KstarToKPiSequence = cms.Sequence(  KstarToKPi  )

BToKstarMuMuSequence = cms.Sequence(
    (muonPairsForKstarMuMu *BToKstarMuMu )
)


BToKstarEESequence = cms.Sequence(
    (electronPairsForKstarEE *BToKstarEE )
)


BToKstarLLSequence = cms.Sequence(
    ( (muonPairsForKstarMuMu *BToKstarMuMu)
     +(electronPairsForKstarEE *BToKstarEE) )   
)


BToKstarLLTables = cms.Sequence( BToKstarEETable + BToKstarMuMuTable )

