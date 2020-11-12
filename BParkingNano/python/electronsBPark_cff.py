import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *

from RecoEgamma.EgammaElectronProducers.lowPtGsfElectronIDExtra_cff import lowPtGsfElectronIDExtra
lowPtGsfElectronExtraID = lowPtGsfElectronIDExtra.clone()
lowPtGsfElectronExtraID.electrons = 'regressionForEle:regressedLowPtElectrons'
lowPtGsfElectronExtraID.rho = 'fixedGridRhoFastjetAll'

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
  #srcMiniAOD              = cms.InputTag('slimmedElectrons',processName=cms.InputTag.skipCurrentProcess()),  
  srcMiniAOD              = cms.InputTag('regressionForEle:regressedElectrons'),
  vertexCollectionMiniAOD = cms.InputTag("offlineSlimmedPrimaryVertices"),
  beamSpotMiniAOD         = cms.InputTag("offlineBeamSpot"),
  conversionsMiniAOD      = cms.InputTag("reducedEgamma:reducedConversions"),
)

electronMVAValueMapProducer = cms.EDProducer(
  'ElectronMVAValueMapProducer',
  # AOD case
  src = cms.InputTag('gedGsfElectrons'),  
  # miniAOD case
  #srcMiniAOD = cms.InputTag('slimmedElectrons',processName=cms.InputTag.skipCurrentProcess()),  
  srcMiniAOD = cms.InputTag('regressionForEle:regressedElectrons'),
    
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

# regression stuff
from RecoEgamma.EgammaTools.regressionModifierNew_cfi import regressionModifier106XUL
from RecoEgamma.EgammaTools.regressionModifierNew_cfi import regressionModifier106XULLP

regressionForEle = cms.EDProducer(
  'ElectronRegresser',
  lowptSrc = cms.InputTag('slimmedLowPtElectrons'),
  pfSrc    = cms.InputTag('slimmedElectrons'),
    lowPtRegressionConfig = cms.PSet(
      modifierName = cms.string('EGRegressionModifierLPV1'),
      rhoTag = cms.string('fixedGridRhoFastjetAll'),
      useClosestToCentreSeedCrysDef = cms.bool(False),
      maxRawEnergyForLowPtEBSigma = cms.double(-1),
      maxRawEnergyForLowPtEESigma = cms.double(1200.),
      eleRegs = cms.PSet(
        ecalOnlyMean = cms.PSet(
            rangeMinLowEt = cms.double(0.2),
            rangeMaxLowEt = cms.double(2.0),
            rangeMinHighEt = cms.double(-1.),
            rangeMaxHighEt = cms.double(3.0),
            forceHighEnergyTrainingIfSaturated = cms.bool(True),
            lowEtHighEtBoundary = cms.double(20.),
            ebLowEtForestName = cms.string("lowPtElectron_eb_ecalOnly_05To20_mean"),
            ebHighEtForestName = cms.string("lowPtElectron_eb_ecalOnly_20To50_mean"),
            eeLowEtForestName = cms.string("lowPtElectron_ee_ecalOnly_05To20_mean"),
            eeHighEtForestName = cms.string("lowPtElectron_ee_ecalOnly_20To50_mean"),
            ),
        ecalOnlySigma = cms.PSet(
            rangeMinLowEt = cms.double(0.0002),
            rangeMaxLowEt = cms.double(0.5),
            rangeMinHighEt = cms.double(0.0002),
            rangeMaxHighEt = cms.double(0.5),
            forceHighEnergyTrainingIfSaturated = cms.bool(True),
            lowEtHighEtBoundary = cms.double(20.),
            ebLowEtForestName = cms.string("lowPtElectron_eb_ecalOnly_05To20_sigma"),
            ebHighEtForestName = cms.string("lowPtElectron_eb_ecalOnly_20To50_sigma"),
            eeLowEtForestName = cms.string("lowPtElectron_ee_ecalOnly_05To20_sigma"),
            eeHighEtForestName = cms.string("lowPtElectron_ee_ecalOnly_20To50_sigma"),
            ),
        epComb = cms.PSet(
            ecalTrkRegressionConfig = cms.PSet(
                rangeMinLowEt = cms.double(0.2),
                rangeMaxLowEt = cms.double(2.0),
                rangeMinHighEt = cms.double(0.2),
                rangeMaxHighEt = cms.double(2.0),
                lowEtHighEtBoundary = cms.double(20.),
                forceHighEnergyTrainingIfSaturated = cms.bool(False),
                ebLowEtForestName = cms.string('lowPtElectron_eb_ecalTrk_05To20_mean'),
                ebHighEtForestName = cms.string('lowPtElectron_eb_ecalTrk_20To50_mean'),
                eeLowEtForestName = cms.string('lowPtElectron_ee_ecalTrk_05To20_mean'),
                eeHighEtForestName = cms.string('lowPtElectron_ee_ecalTrk_20To50_mean'),
                ),
            ecalTrkRegressionUncertConfig = cms.PSet(
                rangeMinLowEt = cms.double(0.0002),
                rangeMaxLowEt = cms.double(0.5),
                rangeMinHighEt = cms.double(0.0002),
                rangeMaxHighEt = cms.double(0.5),
                lowEtHighEtBoundary = cms.double(20.),
                forceHighEnergyTrainingIfSaturated = cms.bool(False),
                ebLowEtForestName = cms.string('lowPtElectron_eb_ecalTrk_05To20_sigma'),
                ebHighEtForestName = cms.string('lowPtElectron_eb_ecalTrk_20To50_sigma'),
                eeLowEtForestName = cms.string('lowPtElectron_ee_ecalTrk_05To20_sigma'),
                eeHighEtForestName = cms.string('lowPtElectron_ee_ecalTrk_20To50_sigma'),
                ),
            maxEcalEnergyForComb=cms.double(200.),
            minEOverPForComb=cms.double(0.025),
            maxEPDiffInSigmaForComb=cms.double(15.),
            maxRelTrkMomErrForComb=cms.double(10.),
            )
        ),
      phoRegs = regressionModifier106XUL.phoRegs.clone()
    ),
    gsfRegressionConfig = cms.PSet(
      modifierName = cms.string('EGRegressionModifierV3'),
      rhoTag = cms.string('fixedGridRhoFastjetAll'),
      useClosestToCentreSeedCrysDef = cms.bool(False),
      maxRawEnergyForLowPtEBSigma = cms.double(-1),
      maxRawEnergyForLowPtEESigma = cms.double(1200.),
      eleRegs = cms.PSet(
        ecalOnlyMean = cms.PSet(
            rangeMinLowEt = cms.double(0.2),
            rangeMaxLowEt = cms.double(2.0),
            rangeMinHighEt = cms.double(-1.),
            rangeMaxHighEt = cms.double(3.0),
            forceHighEnergyTrainingIfSaturated = cms.bool(True),
            lowEtHighEtBoundary = cms.double(999999.),
            ebLowEtForestName = cms.string("gsfElectron_eb_ecalOnly_05To50_mean"),
            ebHighEtForestName = cms.string("gsfElectron_eb_ecalOnly_05To50_mean"),
            eeLowEtForestName = cms.string("gsfElectron_ee_ecalOnly_05To50_mean"),
            eeHighEtForestName = cms.string("gsfElectron_ee_ecalOnly_05To50_mean"),
            ),
        ecalOnlySigma = cms.PSet(
            rangeMinLowEt = cms.double(0.0002),
            rangeMaxLowEt = cms.double(0.5),
            rangeMinHighEt = cms.double(0.0002),
            rangeMaxHighEt = cms.double(0.5),
            forceHighEnergyTrainingIfSaturated = cms.bool(True),
            lowEtHighEtBoundary = cms.double(999999.),
            ebLowEtForestName = cms.string("gsfElectron_eb_ecalOnly_05To50_sigma"),
            ebHighEtForestName = cms.string("gsfElectron_eb_ecalOnly_05To50_sigma"),
            eeLowEtForestName = cms.string("gsfElectron_ee_ecalOnly_05To50_sigma"),
            eeHighEtForestName = cms.string("gsfElectron_ee_ecalOnly_05To50_sigma"),
            ),
        epComb = cms.PSet(
            ecalTrkRegressionConfig = cms.PSet(
                rangeMinLowEt = cms.double(0.2),
                rangeMaxLowEt = cms.double(2.0),
                rangeMinHighEt = cms.double(0.2),
                rangeMaxHighEt = cms.double(2.0),
                lowEtHighEtBoundary = cms.double(999999.),
                forceHighEnergyTrainingIfSaturated = cms.bool(False),
                ebLowEtForestName = cms.string('gsfElectron_eb_ecalTrk_05To50_mean'),
                ebHighEtForestName = cms.string('gsfElectron_eb_ecalTrk_05To50_mean'),
                eeLowEtForestName = cms.string('gsfElectron_ee_ecalTrk_05To50_mean'),
                eeHighEtForestName = cms.string('gsfElectron_ee_ecalTrk_05To50_mean'),
                ),
            ecalTrkRegressionUncertConfig = cms.PSet(
                rangeMinLowEt = cms.double(0.0002),
                rangeMaxLowEt = cms.double(0.5),
                rangeMinHighEt = cms.double(0.0002),
                rangeMaxHighEt = cms.double(0.5),
                lowEtHighEtBoundary = cms.double(999999.),
                forceHighEnergyTrainingIfSaturated = cms.bool(False),
                ebLowEtForestName = cms.string('gsfElectron_eb_ecalTrk_05To50_sigma'),
                ebHighEtForestName = cms.string('gsfElectron_eb_ecalTrk_05To50_sigma'),
                eeLowEtForestName = cms.string('gsfElectron_ee_ecalTrk_05To50_sigma'),
                eeHighEtForestName = cms.string('gsfElectron_ee_ecalTrk_05To50_sigma'),
                ),
            maxEcalEnergyForComb=cms.double(200.),
            minEOverPForComb=cms.double(0.025),
            maxEPDiffInSigmaForComb=cms.double(15.),
            maxRelTrkMomErrForComb=cms.double(10.),
            )
        ),
      phoRegs = regressionModifier106XUL.phoRegs.clone()
    )

)

#Everything can be done here, in one loop and save time :)
electronsForAnalysis = cms.EDProducer(
  'ElectronMerger',
  trgMuon = cms.InputTag('muonTrgSelector:trgMuons'),
  lowptSrc = cms.InputTag('regressionForEle:regressedLowPtElectrons'),
  pfSrc    = cms.InputTag('regressionForEle:regressedElectrons'),
  ptbiasedSeeding = cms.InputTag("lowPtGsfElectronSeedValueMaps","ptbiased","RECO"),
  unbiasedSeeding = cms.InputTag("lowPtGsfElectronSeedValueMaps","unbiased","RECO"),
  mvaId = cms.InputTag("lowPtGsfElectronExtraID"),
  pfmvaId = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2BParkRetrainRawValues"),
  vertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
  ## cleaning wrt trigger muon [-1 == no cut]
  drForCleaning_wrtTrgMuon = cms.double(0.03),
  dzForCleaning_wrtTrgMuon = cms.double(1.),
  ## cleaning between pfEle and lowPtGsf
  drForCleaning = cms.double(0.03),
  dzForCleaning = cms.double(0.5), ##keep tighter dZ to check overlap of pfEle with lowPt (?)
  ## true = flag and clean; false = only flag
  flagAndclean = cms.bool(False),
  pf_ptMin = cms.double(1.),
  ptMin = cms.double(0.5),
  etaMax = cms.double(2.5),
  bdtMin = cms.double(-4.5), #this cut can be used to deactivate low pT e if set to >12
  useRegressionModeForP4 = cms.bool(True),
  useGsfModeForP4 = cms.bool(False),
  sortOutputCollections = cms.bool(True),
  saveLowPtE = cms.bool(True),
    # conversions
    conversions = cms.InputTag('gsfTracksOpenConversions:gsfTracksOpenConversions'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    addUserVarsExtra = cms.bool(True),
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
        convOpen = Var("userInt('convOpen')",bool,doc="Matched to a conversion in gsfTracksOpenConversions collection"),
        convLoose = Var("userInt('convLoose')",bool,doc="Matched to a conversion satisfying Loose WP (see code)"),
        convTight = Var("userInt('convTight')",bool,doc="Matched to a conversion satisfying Tight WP (see code)"),
        convLead = Var("userInt('convLead')",bool,doc="Matched to leading track from conversion"),
        convTrail = Var("userInt('convTrail')",bool,doc="Matched to trailing track from conversion"),
        convExtra = Var("userInt('convExtra')",bool,doc="Flag to indicate if all conversion variables are stored"),
        )
)

if electronsForAnalysis.addUserVarsExtra : 
    electronBParkTable.variables = cms.PSet(
        electronBParkTable.variables,
        convValid = Var("userInt('convValid')",bool,doc="Valid conversion"),
        convChi2Prob = Var("userFloat('convChi2Prob')",float,doc="Reduced chi2 for conversion vertex fit"),
        convQualityHighPurity = Var("userInt('convQualityHighPurity')",bool,doc="'High purity' quality flag for conversion"),
        convQualityHighEff = Var("userInt('convQualityHighEff')",bool,doc="'High efficiency' quality flag for conversion"),
        convTracksN = Var("userInt('convTracksN')",int,doc="Number of tracks associated with conversion"),
        convMinTrkPt = Var("userFloat('convMinTrkPt')",float,doc="Minimum pT found for tracks associated with conversion"),
        convLeadIdx = Var("userInt('convLeadIdx')",int,doc="Index of leading track"),
        convTrailIdx = Var("userInt('convTrailIdx')",int,doc="Index of trailing track"),
        convLxy = Var("userFloat('convLxy')",float,doc="Transverse position of conversion vertex"),
        convVtxRadius = Var("userFloat('convVtxRadius')",float,doc="Radius of conversion vertex"),
        convMass = Var("userFloat('convMass')",float,doc="Invariant mass from conversion pair"),
        convMassFromPin = Var("userFloat('convMassFromPin')",float,doc="Invariant mass from inner momeuntum of conversion pair"),
        convMassBeforeFit = Var("userFloat('convMassBeforeFit')",float,doc="Invariant mass from conversion pair before fit"),
        convMassAfterFit = Var("userFloat('convMassAfterFit')",float,doc="Invariant mass from conversion pair after fit"),
        convLeadNHitsBeforeVtx = Var("userInt('convLeadNHitsBeforeVtx')",int,doc="Number of hits before vertex for lead track"),
        convTrailNHitsBeforeVtx = Var("userInt('convTrailNHitsBeforeVtx')",int,doc="Number of hits before vertex for trail track"),
        convMaxNHitsBeforeVtx = Var("userInt('convMaxNHitsBeforeVtx')",int,doc="Maximum number of hits per track before vertex"),
        convSumNHitsBeforeVtx = Var("userInt('convSumNHitsBeforeVtx')",int,doc="Summed number of hits over tracks before vertex"),
        convDeltaExpectedNHitsInner = Var("userInt('convDeltaExpectedNHitsInner')",int,doc="Delta number of expected hits before vertex"),
        convDeltaCotFromPin = Var("userFloat('convDeltaCotFromPin')",float,doc="Delta cotangent theta from inner momenta"),
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
  regressionForEle
  +lowPtGsfElectronExtraID
  +egmGsfElectronIDSequence
  +electronsForAnalysis
)
electronBParkMC = cms.Sequence(electronsBParkSequence + electronsBParkMCMatchForTable + selectedElectronsMCMatchEmbedded + electronBParkMCTable)
electronBParkTables = cms.Sequence(electronBParkTable)


