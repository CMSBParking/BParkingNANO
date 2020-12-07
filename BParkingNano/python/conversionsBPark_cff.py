import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.BParkingNano.common_cff import *

filteredConversions = cms.EDProducer("ConversionSelector",
    Name = cms.untracked.string('myConvSelection'),
    outFileName = cms.string('ValidationHistos.root'),
    convProducer = cms.string('gsfTracksOpenConversions'),
    conversionCollection = cms.string('gsfTracksOpenConversions'),
    lowPtElectronProducer = cms.string('slimmedLowPtElectrons'),
    lowPtElectronCollection =  cms.string(''),
    unbiasedSeeding = cms.InputTag("lowPtGsfElectronSeedValueMaps","unbiased","RECO"),
    mvaId = cms.InputTag("lowPtGsfElectronLatestID"),
    phoProducer = cms.string('photons'),
    photonCollection = cms.string(''),
    filteredConversions = cms.string(''),
    rhoHandle = cms.string('fixedGridRhoAll'),                                    
    dqmpath = cms.string('EgammaV/ConversionValidator/'),
    Verbosity = cms.untracked.int32(0),
    generalTracksOnly = cms.bool(False),
    gsfTracksOpenOnly = cms.bool(True),
    arbitratedMerged =  cms.bool(False),
    arbitratedEcalSeeded = cms.bool(False),
    ecalalgotracks = cms.bool(False),
    minPtOnTracks = cms.double(0.5),
    highPurity = cms.bool(True),
    minProb = cms.double(0.0005),
    maxHitsBeforeVtx = cms.uint32(1),
    minLxy = cms.double(0.),

    minPhoPtForEffic = cms.double(0.3),#when hardcoded it was 2.5
    maxPhoEtaForEffic = cms.double(2.5),
    maxPhoZForEffic = cms.double(200.),
    maxPhoRForEffic = cms.double(100.),
    minPhoPtForPurity = cms.double(0.),#when hardcoded it was 0.5
    maxPhoEtaForPurity = cms.double(2.5),
    maxPhoZForPurity = cms.double(200.),
    maxPhoRForPurity = cms.double(100.),

#
    minPhoEtCut = cms.double(0.),

#
    useTP =  cms.bool(True),
#

    etBin = cms.int32(100),                                  
    etMax = cms.double(100.0),                                  
    etMin = cms.double(0.0),
#
    etaBin = cms.int32(100),
    etaBin2 = cms.int32(25),
    etaMin = cms.double(-2.5),
    etaMax = cms.double(2.5),
#
    phiBin = cms.int32(100),
    phiMin = cms.double(-3.14),
    phiMax = cms.double(3.14),
#
    resBin = cms.int32(100),
    resMin = cms.double(0.),
    resMax = cms.double(1.1),
#
    eoverpBin =  cms.int32(100),
    eoverpMin =  cms.double(0.),
    eoverpMax =  cms.double(5.),
#                                       
    dEtaTracksBin = cms.int32(100),
    dEtaTracksMin = cms.double(-0.2),
    dEtaTracksMax = cms.double(0.2),
#
    dPhiTracksBin = cms.int32(100),
    dPhiTracksMin = cms.double(-0.5),
    dPhiTracksMax = cms.double(0.5),
#
    dEtaBin = cms.int32(100),
    dEtaMin = cms.double(-0.2),
    dEtaMax = cms.double(0.2),
#
    dPhiBin = cms.int32(100),
    dPhiMin = cms.double(-0.05),
    dPhiMax = cms.double(0.05),
#
    rBin = cms.int32(60), 
    rMin = cms.double(0.),
    rMax = cms.double(120),
#
    zBin = cms.int32(100),
    zMin = cms.double(-220.),
    zMax = cms.double(220),
#
 
    dCotTracksBin = cms.int32(100),                              
    dCotTracksMin = cms.double(-0.12),
    dCotTracksMax = cms.double(0.12),
#                                  
    chi2Min =  cms.double(0.),
    chi2Max =  cms.double(20.),                              
#

    rBinForXray = cms.int32(200),
    rMinForXray = cms.double(0.),
    rMaxForXray = cms.double(80.),                               
    zBinForXray = cms.int32(100),
    zBin2ForXray = cms.int32(560),
    zMinForXray = cms.double(0.),
    zMaxForXray = cms.double(280.),                               

    simTracks = cms.InputTag("g4SimHits")
)


conversionTable  = cms.EDProducer("SimpleCompositeCandidateFlatTableProducer",
    #src = cms.InputTag("gsfTracksOpenConversions","gsfTracksOpenConversions"),
    src = cms.InputTag("filteredConversions"),
    cut = cms.string(""), #we should not filter on cross linked collections
    name = cms.string("Conversion"),
    doc  = cms.string("Our packedCandidate table"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the muons
    variables = cms.PSet(CandVars,
                         rho = ufloat('rho'),
                         vtxChi2Prob = ufloat('vtxChi2Prob'),
                         vtxR   = ufloat('vtxR'),
                         vtxX   = ufloat('vtxX'),
                         vtxY   = ufloat('vtxY'),
                         vtxZ   = ufloat('vtxZ'),
                         dCotFromPin = ufloat('dCotFromPin'),
                         dMinAppTracks = ufloat('dMinAppTracks'),
                         dPhiTracksAtVtx = ufloat('dPhiTracksAtVtx'),
# 
                         tkLead_pt = ufloat('tkLead_pt'),
                         tkTrail_pt = ufloat('tkTrail_pt'),
                         tkLead_eta = ufloat('tkLead_eta'),
                         tkTrail_eta = ufloat('tkTrail_eta'),
                         tkLead_phi = ufloat('tkLead_phi'),
                         tkTrail_phi = ufloat('tkTrail_phi'),
                         tkLead_nHits = ufloat('tkLead_nHits'),
                         tkTrail_nHits = ufloat('tkTrail_nHits'),
                         tkLead_isMatchedToLowPtEle = uint('tkLead_isMatchedToLowPtEle'),
                         tkTrail_isMatchedToLowPtEle = uint('tkTrail_isMatchedToLowPtEle'),
#
                         eleLead_pt = ufloat('eleLead_pt'),
                         eleTrail_pt = ufloat('eleTrail_pt'),
                         eleLead_ch = ufloat('eleLead_ch'),
                         eleTrail_ch = ufloat('eleTrail_ch'),
                         eleLead_eta = ufloat('eleLead_eta'),
                         eleTrail_eta = ufloat('eleTrail_eta'),
                         eleLead_phi = ufloat('eleLead_phi'),
                         eleTrail_phi = ufloat('eleTrail_phi'),
#
                         eleLead_scEnergy = ufloat('eleLead_scEnergy'),
                         eleTrail_scEnergy = ufloat('eleTrail_scEnergy'),
                         eleLead_scEta = ufloat('eleLead_scEta'),
                         eleTrail_scEta = ufloat('eleTrail_scEta'),
                         eleLead_scEtaWidth = ufloat('eleLead_scEtaWidth'),
                         eleTrail_scEtaWidth = ufloat('eleTrail_scEtaWidth'),
                         eleLead_scPhiWidth = ufloat('eleLead_scPhiWidth'),
                         eleTrail_scPhiWidth = ufloat('eleTrail_scPhiWidth'),
                         eleLead_cluSize = ufloat('eleLead_cluSize'),
                         eleTrail_cluSize = ufloat('eleTrail_cluSize'),
                         eleLead_full5x5_HoverE = ufloat('eleLead_full5x5_HoverE'),
                         eleTrail_full5x5_HoverE = ufloat('eleTrail_full5x5_HoverE'),
                         eleLead_full5x5_r9 = ufloat('eleLead_full5x5_r9'),
                         eleTrail_full5x5_r9 = ufloat('eleTrail_full5x5_r9'),
#
                         eleLead_seed_dEta = ufloat('eleLead_seed_dEta'),
                         eleTrail_seed_dEta = ufloat('eleTrail_seed_dEta'),
                         eleLead_eclu_EOverP =ufloat('eleLead_eclu_EOverP'),
                         eleTrail_eclu_EOverP =ufloat('eleTrail_eclu_EOverP'),
               
                         eleLead_scEOverP = ufloat('eleLead_scEOverP'),
                         eleTrail_scEOverP = ufloat('eleTrail_scEOverP'),
                         eleLead_sc_dEta = ufloat('eleLead_sc_dEta'),
                         eleTrail_sc_dEta = ufloat('eleTrail_sc_dEta'),
                         eleLead_sc_dPhi = ufloat('eleLead_sc_dPhi'),
                         eleTrail_sc_dPhi = ufloat('eleTrail_sc_dPhi'),

                         eleLead_fBrem = ufloat('eleLead_fBrem'),
                         eleTrail_fBrem = ufloat('eleTrail_fBrem'),
                         eleLead_shFracHits = ufloat('eleLead_shFracHits'),
                         eleTrail_shFracHits = ufloat('eleTrail_shFracHits'),

                         eleLead_pMode = ufloat('eleLead_pMode'),
                         eleTrail_pMode = ufloat('eleTrail_pMode'),
                         eleLead_chi2 = ufloat('eleLead_chi2'),
                         eleTrail_chi2 = ufloat('eleTrail_chi2'),
                         eleLead_nHits = ufloat('eleLead_nHits'),
                         eleTrail_nHits = ufloat('eleTrail_nHits'),
                         eleLead_dR = ufloat('eleLead_dR'),
                         eleTrail_dR = ufloat('eleTrail_dR'),
#
                         eleLead_unbiasedSeedBDT = ufloat('eleLead_unbiasedSeedBDT'),
                         eleTrail_unbiasedSeedBDT = ufloat('eleTrail_unbiasedSeedBDT'),
                         eleLead_mvaId = ufloat('eleLead_mvaId'),
                         eleTrail_mvaId = ufloat('eleTrail_mvaId'),






        )
)

from Configuration.Eras.Modifier_fastSim_cff import fastSim
if fastSim.isChosen():
    tkConversionValidation.simTracks = cms.InputTag("famosSimHits")

#conversionBParkSequence = cms.Sequence(filteredConversions+conversionTable)
conversionBParkSequence = cms.Sequence(filteredConversions)
