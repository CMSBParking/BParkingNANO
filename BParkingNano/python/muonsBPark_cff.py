import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *



muonTrgSelector = cms.EDProducer("MuonTriggerSelector",
                                 muonCollection = cms.InputTag("slimmedMuons"), #same collection as in NanoAOD                                                           
                                 bits = cms.InputTag("TriggerResults","","HLT"),
                                 prescales = cms.InputTag("patTrigger"),
                                 objects = cms.InputTag("slimmedPatTrigger"),
                                 vertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                 
                                 ##for the output trigger matched collection
                                 maxdR_matching = cms.double(0.1),
                                 
                                 ## for the output selected collection (tag + all compatible in dZ)
                                 dzForCleaning_wrtTrgMuon = cms.double(1.),

                                 ptMin = cms.double(0.5),
                                 absEtaMax = cms.double(2.4),
                                 # keeps only muons with at soft Quality flag
                                 softMuonsOnly = cms.bool(False)
                             )

countTrgMuons = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("muonTrgSelector", "trgMuons")
)


muonBParkTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
    src = cms.InputTag("muonTrgSelector:SelectedMuons"),
    cut = cms.string(""), #we should not filter on cross linked collections
    name = cms.string("Muon"),
    doc  = cms.string("slimmedMuons for BPark after basic selection"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the muons
    variables = cms.PSet(CandVars,
        ptErr   = Var("bestTrack().ptError()", float, doc = "ptError of the muon track", precision=6),
        dz = Var("dB('PVDZ')",float,doc="dz (with sign) wrt first PV, in cm",precision=10),
        dzErr = Var("abs(edB('PVDZ'))",float,doc="dz uncertainty, in cm",precision=6),
        dxy = Var("dB('PV2D')",float,doc="dxy (with sign) wrt first PV, in cm",precision=10),
        dxyErr = Var("edB('PV2D')",float,doc="dxy uncertainty, in cm",precision=6),
        vx = Var("vx()",float,doc="x coordinate of vertex position, in cm",precision=6),
        vy = Var("vy()",float,doc="y coordinate of vertex position, in cm",precision=6),
        vz = Var("vz()",float,doc="z coordinate of vertex position, in cm",precision=6),
        ip3d = Var("abs(dB('PV3D'))",float,doc="3D impact parameter wrt first PV, in cm",precision=10),
        sip3d = Var("abs(dB('PV3D')/edB('PV3D'))",float,doc="3D impact parameter significance wrt first PV",precision=10),
#        segmentComp   = Var("segmentCompatibility()", float, doc = "muon segment compatibility", precision=14), # keep higher precision since people have cuts with 3 digits on this
#        nStations = Var("numberOfMatchedStations", int, doc = "number of matched stations with default arbitration (segment & track)"),
        #nTrackerLayers = Var("innerTrack().hitPattern().trackerLayersWithMeasurement()", int, doc = "number of layers in the tracker"),
#        pfRelIso03_chg = Var("pfIsolationR03().sumChargedHadronPt/pt",float,doc="PF relative isolation dR=0.3, charged component"),
        pfRelIso03_all = Var("(pfIsolationR03().sumChargedHadronPt + max(pfIsolationR03().sumNeutralHadronEt + pfIsolationR03().sumPhotonEt - pfIsolationR03().sumPUPt/2,0.0))/pt",float,doc="PF relative isolation dR=0.3, total (deltaBeta corrections)"),
        pfRelIso04_all = Var("(pfIsolationR04().sumChargedHadronPt + max(pfIsolationR04().sumNeutralHadronEt + pfIsolationR04().sumPhotonEt - pfIsolationR04().sumPUPt/2,0.0))/pt",float,doc="PF relative isolation dR=0.4, total (deltaBeta corrections)"),
#        tightCharge = Var("?(muonBestTrack().ptError()/muonBestTrack().pt() < 0.2)?2:0",int,doc="Tight charge criterion using pterr/pt of muonBestTrack (0:fail, 2:pass)"),
        isPFcand = Var("isPFMuon",bool,doc="muon is PF candidate"),
        isGlobal = Var("isGlobalMuon",bool,doc="muon is global muon"),
        isTracker = Var("isTrackerMuon",bool,doc="muon is tracker muon"),
        mediumId = Var("passed('CutBasedIdMedium')",bool,doc="cut-based ID, medium WP"),
#        mediumPromptId = Var("passed('CutBasedIdMediumPrompt')",bool,doc="cut-based ID, medium prompt WP"),
        tightId = Var("passed('CutBasedIdTight')",bool,doc="cut-based ID, tight WP"),
        softId = Var("passed('SoftCutBasedId')",bool,doc="soft cut-based ID"),
#        softMvaId = Var("passed('SoftMvaId')",bool,doc="soft MVA ID"),
#        highPtId = Var("?passed('CutBasedIdGlobalHighPt')?2:passed('CutBasedIdTrkHighPt')","uint8",doc="high-pT cut-based ID (1 = tracker high pT, 2 = global high pT, which includes tracker high pT)"),
        pfIsoId = Var("passed('PFIsoVeryLoose')+passed('PFIsoLoose')+passed('PFIsoMedium')+passed('PFIsoTight')+passed('PFIsoVeryTight')+passed('PFIsoVeryVeryTight')","uint8",doc="PFIso ID from miniAOD selector (1=PFIsoVeryLoose, 2=PFIsoLoose, 3=PFIsoMedium, 4=PFIsoTight, 5=PFIsoVeryTight, 6=PFIsoVeryVeryTight)"),
        tkIsoId = Var("?passed('TkIsoTight')?2:passed('TkIsoLoose')","uint8",doc="TkIso ID (1=TkIsoLoose, 2=TkIsoTight)"),
#        mvaId = Var("passed('MvaLoose')+passed('MvaMedium')+passed('MvaTight')","uint8",doc="Mva ID from miniAOD selector (1=MvaLoose, 2=MvaMedium, 3=MvaTight)"),
#        miniIsoId = Var("passed('MiniIsoLoose')+passed('MiniIsoMedium')+passed('MiniIsoTight')+passed('MiniIsoVeryTight')","uint8",doc="MiniIso ID from miniAOD selector (1=MiniIsoLoose, 2=MiniIsoMedium, 3=MiniIsoTight, 4=MiniIsoVeryTight)"),
#        multiIsoId = Var("?passed('MultiIsoMedium')?2:passed('MultiIsoLoose')","uint8",doc="MultiIsoId from miniAOD selector (1=MultiIsoLoose, 2=MultiIsoMedium)"),
        triggerIdLoose = Var("passed('TriggerIdLoose')",bool,doc="TriggerIdLoose ID"),
#        inTimeMuon = Var("passed('InTimeMuon')",bool,doc="inTimeMuon ID"),
        isTriggering = Var("userInt('isTriggering')", int,doc="flag the reco muon is also triggering")
    ),
)


muonsBParkMCMatchForTable = cms.EDProducer("MCMatcher",            # cut on deltaR, deltaPt/Pt; pick best by deltaR
    src         = muonBParkTable.src,                         # final reco collection
    matched     = cms.InputTag("finalGenParticlesBPark"),     # final mc-truth particle collection
    mcPdgId     = cms.vint32(13),                             # one or more PDG ID (13 = mu); absolute values (see below)
    checkCharge = cms.bool(False),                            # True = require RECO and MC objects to have the same charge
    mcStatus    = cms.vint32(1),                              # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
    maxDeltaR   = cms.double(0.03),                           # Minimum deltaR for the match
    maxDPtRel   = cms.double(0.5),                            # Minimum deltaPt/Pt for the match
    resolveAmbiguities    = cms.bool(True),                   # Forbid two RECO objects to match to the same GEN object
    resolveByMatchQuality = cms.bool(True),                   # False = just match input in order; True = pick lowest deltaR pair first
)

muonBParkMCTable = cms.EDProducer("CandMCMatchTableProducerBPark",
    src     = muonBParkTable.src,
    mcMap   = cms.InputTag("muonsBParkMCMatchForTable"),
    objName = muonBParkTable.name,
    objType = muonBParkTable.name, 
    branchName = cms.string("genPart"),
    docString = cms.string("MC matching to status==1 muons"),
)

selectedMuonsMCMatchEmbedded = cms.EDProducer(
    'MuonMatchEmbedder',
    src = cms.InputTag('muonTrgSelector:SelectedMuons'),
    matching = cms.InputTag('muonsBParkMCMatchForTable')
)

muonTriggerMatchedTable = muonBParkTable.clone(
    src = cms.InputTag("muonTrgSelector:trgMuons"),
    name = cms.string("TriggerMuon"),
    doc  = cms.string("reco muon matched to triggering muon"),
    variables = cms.PSet(CandVars,
        vx = Var("vx()",float,doc="x coordinate of vertex position, in cm",precision=6),
        vy = Var("vy()",float,doc="y coordinate of vertex position, in cm",precision=6),
        vz = Var("vz()",float,doc="z coordinate of vertex position, in cm",precision=6),
        trgMuonIndex = Var("userInt('trgMuonIndex')", int,doc="index in trigger muon collection")
   )
)

muonBParkSequence = cms.Sequence(muonTrgSelector * countTrgMuons)
muonBParkMC = cms.Sequence(muonBParkSequence + muonsBParkMCMatchForTable + selectedMuonsMCMatchEmbedded + muonBParkMCTable)
muonBParkTables = cms.Sequence(muonBParkTable)
muonTriggerMatchedTables = cms.Sequence(muonTriggerMatchedTable)
