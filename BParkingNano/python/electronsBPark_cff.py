import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *

from RecoEgamma.EgammaElectronProducers.lowPtGsfElectronID_cff import lowPtGsfElectronID
lowPtGsfElectronLatestID = lowPtGsfElectronID.clone()
lowPtGsfElectronLatestID.electrons = 'slimmedLowPtElectrons'
lowPtGsfElectronLatestID.rho = 'fixedGridRhoFastjetAll'

##essentially the commented out can be inside the same loop... no need to have a more loops in an "expensive" object
'''lowptElectronsWithSeed = cms.EDProducer(
  'PATLowPtElectronSeedingEmbedder',
  src = cms.InputTag('slimmedLowPtElectrons'),
  ptbiasedSeeding = cms.InputTag("lowPtGsfElectronSeedValueMaps","ptbiased","RECO"),
  unbiasedSeeding = cms.InputTag("lowPtGsfElectronSeedValueMaps","unbiased","RECO"),
    minBdtUnbiased = cms.double(0.5)
)
lowptElectronsForAnalysis = cms.EDFilter(
  'PATElectronSelector',
  src = cms.InputTag("lowptElectronsWithSeed"),
  ## need to add cut on BDT and ID when available 
  ## pT > 0.5 to accomodate l1 and l2
  cut = cms.string('pt > 0.5 && eta > -2.4 && eta < 2.4'),
  )
pfElectronsForAnalysis = cms.EDFilter(
  'PATElectronSelector',
  src = cms.InputTag("slimmedElectrons"),
  ## can add cut on ID.
  ## pT > 2 since almost nothing below anyway
  cut = cms.string("pt > 2 && eta > -2.4 && eta < 2.4"),
  )
'''




#Everything can be done here, in one loop and save time :)
electronsForAnalysis = cms.EDProducer(
  'ElectronMerger',
  trgMuon = cms.InputTag('muonTrgSelector:trgMuons'),
  lowptSrc = cms.InputTag('slimmedLowPtElectrons'),
  pfSrc    = cms.InputTag('slimmedElectrons'),
  ptbiasedSeeding = cms.InputTag("lowPtGsfElectronSeedValueMaps","ptbiased","RECO"),
  unbiasedSeeding = cms.InputTag("lowPtGsfElectronSeedValueMaps","unbiased","RECO"),
  mvaId = cms.InputTag("lowPtGsfElectronLatestID"),
  vertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
  ## cleaning wrt trigger muon [-1 == no cut]
  drForCleaning_wrtTrgMuon = cms.double(0.03),
  dzForCleaning_wrtTrgMuon = cms.double(1.),
  ## cleaning between pfEle and lowPtGsf
  drForCleaning = cms.double(0.03),
  dzForCleaning = cms.double(0.7), ##keep tighter dZ to check overlap of pfEle with lowPt (?)
  ## true = flag and clean; false = only flag
  flagAndclean = cms.bool(False),
  pf_ptMin = cms.double(1.),
  ptMin = cms.double(0.5),
  etaMax = cms.double(2.5),
  bdtMin = cms.double(-4), #this cut can be used to deactivate low pT e if set to >12
  useGsfModeForP4 = cms.bool(True),
  sortOutputCollections = cms.bool(True)
)

electronBParkTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
 src = cms.InputTag("electronsForAnalysis:SelectedElectrons"),
 cut = cms.string(""),
    name= cms.string("Electron"),
    doc = cms.string("slimmedElectrons for BPark after basic selection"),
    singleton = cms.bool(False), 
    extension = cms.bool(False),                                                
    variables = cms.PSet(P4Vars,
        pdgId  = Var("pdgId", int, doc="PDG code assigned by the event reconstruction (not by MC truth)"),
        charge = Var("userFloat('chargeMode')", int, doc="electric charge from pfEle or chargeMode for lowPtGsf"),
        dz = Var("dB('PVDZ')",float,doc="dz (with sign) wrt first PV, in cm",precision=10),
        dzErr = Var("abs(edB('PVDZ'))",float,doc="dz uncertainty, in cm",precision=6),
        dxy = Var("dB('PV2D')",float,doc="dxy (with sign) wrt first PV, in cm",precision=10),
        dxyErr = Var("edB('PV2D')",float,doc="dxy uncertainty, in cm",precision=6),
        vx = Var("vx()",float,doc="x coordinate of vertex position, in cm",precision=6),
        vy = Var("vy()",float,doc="y coordinate of vertex position, in cm",precision=6),
        vz = Var("vz()",float,doc="z coordinate of vertex position, in cm",precision=6),
        ip3d = Var("abs(dB('PV3D'))",float,doc="3D impact parameter wrt first PV, in cm",precision=10),
        sip3d = Var("abs(dB('PV3D')/edB('PV3D'))",float,doc="3D impact parameter significance wrt first PV, in cm",precision=10),
        deltaEtaSC = Var("superCluster().eta()-eta()",float,doc="delta eta (SC,ele) with sign",precision=10),
        r9 = Var("full5x5_r9()",float,doc="R9 of the supercluster, calculated with full 5x5 region",precision=10),
        sieie = Var("full5x5_sigmaIetaIeta()",float,doc="sigma_IetaIeta of the supercluster, calculated with full 5x5 region",precision=10),
        hoe = Var("hadronicOverEm()",float,doc="H over E",precision=8),
        tightCharge = Var("isGsfCtfScPixChargeConsistent() + isGsfScPixChargeConsistent()",int,doc="Tight charge criteria (0:none, 1:isGsfScPixChargeConsistent, 2:isGsfCtfScPixChargeConsistent)"),
        convVeto = Var("passConversionVeto()",bool,doc="pass conversion veto"),
        lostHits = Var("gsfTrack.hitPattern.numberOfLostHits('MISSING_INNER_HITS')","uint8",doc="number of missing inner hits"),
        isPF = Var("userInt('isPF')",bool,doc="electron is PF candidate"),
        isLowPt = Var("userInt('isLowPt')",bool,doc="electron is LowPt candidate"),
        ptBiased = Var("userFloat('ptBiased')",float,doc="ptBiased from seed BDT 20 for pfEle"), 
        unBiased = Var("userFloat('unBiased')",float,doc="unBiased from seed BDT 20 for pfEle"), 
        mvaId = Var("userFloat('mvaId')",float,doc="MVA ID for low pT, 20 for pfEle"),
        fBrem = Var("fbrem()",float,doc="brem fraction from the gsf fit",precision=8),
        isPFoverlap = Var("userInt('isPFoverlap')",bool,doc="flag lowPt ele overlapping with pf in selected_pf_collection",precision=8),
        )
)



electronsBParkMCMatchForTable = cms.EDProducer("MCMatcher",  # cut on deltaR, deltaPt/Pt; pick best by deltaR
    src         = electronBParkTable.src,                 # final reco collection
    matched     = cms.InputTag("finalGenParticlesBPark"), # final mc-truth particle collection
    mcPdgId     = cms.vint32(11,22),                 # one or more PDG ID (11 = el, 22 = pho); absolute values (see below)
    checkCharge = cms.bool(False),              # True = require RECO and MC objects to have the same charge  
    mcStatus    = cms.vint32(1),                # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
    maxDeltaR   = cms.double(0.03),             # Minimum deltaR for the match
    maxDPtRel   = cms.double(0.5),              # Minimum deltaPt/Pt for the match
    resolveAmbiguities    = cms.bool(False),    # Forbid two RECO objects to match to the same GEN object
    resolveByMatchQuality = cms.bool(True),    # False = just match input in order; True = pick lowest deltaR pair first
    
)

selectedElectronsMCMatchEmbedded = cms.EDProducer(
    'ElectronMatchEmbedder',
    src = electronBParkTable.src,
    matching = cms.InputTag('electronsBParkMCMatchForTable')
)

electronBParkMCTable = cms.EDProducer("CandMCMatchTableProducerBPark",
    src     = electronBParkTable.src,
    mcMap   = cms.InputTag("electronsBParkMCMatchForTable"),
    objName = electronBParkTable.name,
    objType = electronBParkTable.name,
    branchName = cms.string("genPart"),
    docString = cms.string("MC matching to status==1 electrons or photons"),
)


electronsBParkSequence = cms.Sequence(
  lowPtGsfElectronLatestID
  +electronsForAnalysis
)
electronBParkMC = cms.Sequence(electronsBParkSequence + electronsBParkMCMatchForTable + selectedElectronsMCMatchEmbedded + electronBParkMCTable)
electronBParkTables = cms.Sequence(electronBParkTable)


