import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *

from RecoEgamma.EgammaElectronProducers.lowPtGsfElectronID_cff import lowPtGsfElectronID
lowPtGsfElectronLatestID = lowPtGsfElectronID.clone()
lowPtGsfElectronLatestID.electrons = 'slimmedLowPtElectrons'
lowPtGsfElectronLatestID.rho = 'fixedGridRhoFastjetAll'


mvaConfigsForEleProducer = cms.VPSet( )
# Import and add all desired MVAs
from PhysicsTools.BParkingNano.mvaElectronID_BParkRetrain_cff \
    import mvaEleID_BParkRetrain_producer_config
mvaConfigsForEleProducer.append( mvaEleID_BParkRetrain_producer_config )


# The producer to compute the MVA input variables which are not accessible with the cut parser
electronMVAVariableHelper = cms.EDProducer('GsfElectronMVAVariableHelper',
  # The module automatically detects AOD vs miniAOD, so we configure both
  # AOD case
  src = cms.InputTag('gedGsfElectrons'),
  vertexCollection = cms.InputTag("offlinePrimaryVertices"),
  beamSpot         = cms.InputTag("offlineBeamSpot"),
conversions = cms.InputTag("allConversions"),
  # miniAOD case
  srcMiniAOD              = cms.InputTag('slimmedElectrons',processName=cms.InputTag.skipCurrentProcess()),
  vertexCollectionMiniAOD = cms.InputTag("offlineSlimmedPrimaryVertices"),
  beamSpotMiniAOD         = cms.InputTag("offlineBeamSpot"),
  conversionsMiniAOD      = cms.InputTag("reducedEgamma:reducedConversions"),
)

electronMVAValueMapProducer = cms.EDProducer(
  'ElectronMVAValueMapProducer',
  # AOD case
  src = cms.InputTag('gedGsfElectrons'),  
  # miniAOD case
  srcMiniAOD = cms.InputTag('slimmedElectrons',processName=cms.InputTag.skipCurrentProcess()),
    
  # MVA configurations
  mvaConfigurations = mvaConfigsForEleProducer
)

egmGsfElectronIDs = cms.EDProducer(
    "VersionedGsfElectronIdProducer",
    physicsObjectSrc = cms.InputTag('gedGsfElectrons'),
    physicsObjectIDs = cms.VPSet( )
)

egmGsfElectronIDTask = cms.Task(
    electronMVAVariableHelper,
    electronMVAValueMapProducer,
    egmGsfElectronIDs,
)
egmGsfElectronIDSequence = cms.Sequence(egmGsfElectronIDTask)


SaveLowPtE=True #Flag to store the LowpT e
OnlyPFeCompatible=False #flag to take care compatibility when only PF e is availiable (eg 2017 Run)


if OnlyPFeCompatible:
   lowpt_input='slimmedElectrons'
   biasedseed_input='electronMVAValueMapProducer:ElectronMVAEstimatorRun2BParkRetrainRawValues'
   unbiasedseed_input='electronMVAValueMapProducer:ElectronMVAEstimatorRun2BParkRetrainRawValues'
   id_input='electronMVAValueMapProducer:ElectronMVAEstimatorRun2BParkRetrainRawValues'
   SaveLowPtE=False
else:
   lowpt_input='slimmedLowPtElectrons'
   biasedseed_input='lowPtGsfElectronSeedValueMaps:ptbiased'
   unbiasedseed_input='lowPtGsfElectronSeedValueMaps:unbiased'
   id_input='lowPtGsfElectronLatestID'

#Everything can be done here, in one loop and save time :)
electronsForAnalysis = cms.EDProducer(
  'ElectronMerger',
  trgMuon = cms.InputTag('muonTrgSelector:trgMuons'),
  lowptSrc = cms.InputTag(lowpt_input),
  pfSrc    = cms.InputTag('slimmedElectrons'),
  ptbiasedSeeding = cms.InputTag(biasedseed_input),
  unbiasedSeeding = cms.InputTag(unbiasedseed_input),
  mvaId = cms.InputTag(id_input),
  pfmvaId = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2BParkRetrainRawValues"),
  vertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
  ## cleaning wrt trigger muon [-1 == no cut]
  drForCleaning_wrtTrgMuon = cms.double(0.03), #was positive
  dzForCleaning_wrtTrgMuon = cms.double(1.),  #was positive 1
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
  sortOutputCollections = cms.bool(True),
  saveLowPtE = cms.bool(SaveLowPtE) #skips low pT electrons

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
#        deltaEtaSC = Var("superCluster().eta()-eta()",float,doc="delta eta (SC,ele) with sign",precision=10),
#        r9 = Var("full5x5_r9()",float,doc="R9 of the supercluster, calculated with full 5x5 region",precision=10),
#        sieie = Var("full5x5_sigmaIetaIeta()",float,doc="sigma_IetaIeta of the supercluster, calculated with full 5x5 region",precision=10),
#        hoe = Var("hadronicOverEm()",float,doc="H over E",precision=8),
#        tightCharge = Var("isGsfCtfScPixChargeConsistent() + isGsfScPixChargeConsistent()",int,doc="Tight charge criteria (0:none, 1:isGsfScPixChargeConsistent, 2:isGsfCtfScPixChargeConsistent)"),
        convVeto = Var("passConversionVeto()",bool,doc="pass conversion veto"),
#        lostHits = Var("gsfTrack.hitPattern.numberOfLostHits('MISSING_INNER_HITS')","uint8",doc="number of missing inner hits"),
        pfRelIso = Var("(pfIsolationVariables().sumChargedHadronPt+max(0.0,pfIsolationVariables().sumNeutralHadronEt+pfIsolationVariables().sumPhotonEt-0.5*pfIsolationVariables().sumPUPt))/pt",float,doc="PF relative isolation dR=0.3, total (deltaBeta corrections)"),
        trkRelIso = Var("trackIso/pt",float,doc="PF relative isolation dR=0.3, total (deltaBeta corrections)"),
        isPF = Var("userInt('isPF')",bool,doc="electron is PF candidate"),
        isLowPt = Var("userInt('isLowPt')",bool,doc="electron is LowPt candidate"),
        ptBiased = Var("userFloat('ptBiased')",float,doc="ptBiased from seed BDT 20 for pfEle"), 
        unBiased = Var("userFloat('unBiased')",float,doc="unBiased from seed BDT 20 for pfEle"), 
        mvaId = Var("userFloat('mvaId')",float,doc="MVA ID for low pT, 20 for pfEle"),
        pfmvaId = Var("userFloat('pfmvaId')",float,doc="MVA ID for pfEle, 20 for low pT"),
        fBrem = Var("fbrem()",float,doc="brem fraction from the gsf fit",precision=12),
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

if SaveLowPtE:
  electronsBParkSequence = cms.Sequence(
     lowPtGsfElectronLatestID
     +egmGsfElectronIDSequence
     +electronsForAnalysis
   )
else:
  electronsBParkSequence = cms.Sequence(
     egmGsfElectronIDSequence
     +electronsForAnalysis
   )

electronBParkMC = cms.Sequence(electronsBParkSequence + electronsBParkMCMatchForTable + selectedElectronsMCMatchEmbedded + electronBParkMCTable)
electronBParkTables = cms.Sequence(electronBParkTable)


