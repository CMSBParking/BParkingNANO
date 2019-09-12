import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *

muonsTracksBPark = cms.EDProducer('MuonTrackMerger',
                                  muons  = cms.InputTag("muonTrgSelector:SelectedMuons"),
                                  muonTransientTracks = cms.InputTag("muonTrgSelector:SelectedTransientMuons"),
                                  tracks = cms.InputTag("tracksBPark:SelectedTracks"),
                                  trackTransientTracks = cms.InputTag("tracksBPark:SelectedTransientTracks"),
                                  muonSelection = cms.string(''),
                                  ## pt vs eta range from AN 2008/097 Table2
                                  trackSelection = cms.string('( userInt("isPacked") == 1 && userInt("nValidHits") > 15 && userInt("nValidHits") < 20 && ( (pt() < 4.6 && eta() > -1.2 && eta() < 1.2 ) || (pt() < 3.6 && ( (eta() > -1.5 && eta() < -1.2) || (eta() > 1.2 && eta() < 1.5) ) ) || (pt() < 1.2 && ( (eta() > -2.4 && eta() < 1.5 ) || (eta() > 1.5 && eta() < 2.4) ) ) ) ) && userInt("isMatchedToMuon") == 0 && userInt("isMatchedToEle") == 0'),
                                  sortOutputCollections = cms.bool(True)
                                  )


muonsTracksBParkTable = cms.EDProducer(
    "SimpleCompositeCandidateFlatTableProducer",
    src = cms.InputTag("muonsTracksBPark:SelectedMuonsTracks"),
    cut = cms.string(""),
    name = cms.string("muonTracks"),
    doc  = cms.string("pfMuon and muon-like tracks in the EB acceptance region"),
    singleton = cms.bool(False),
    extension = cms.bool(False),
    variables = cms.PSet(
         CandVars,
        vx = Var("vx()", float, doc="x coordinate of vertex position, in cm", precision=10),
        vy = Var("vy()", float, doc="y coordinate of vertex position, in cm", precision=10),
        vz = Var("vz()", float, doc="z coordinate of vertex position, in cm", precision=10),
        isTrack = Var("userInt('isTrack')",int,doc="muon-like track", precision=10),
        isMuon = Var("userInt('isMuon')",int,doc="genuine muon from MuonSelected collection", precision=10),
        dz = Var("userFloat('dz')",float,doc="dz (with sign) wrt first PV, in cm", precision=10),
        dxy = Var("userFloat('dxy')",float,doc="dxy (with sign) wrt first PV, in cm", precision=10),
        dzS = Var("userFloat('dzS')", float, doc="dz/err (with sign) wrt first PV, in cm", precision=10),
        dxyS = Var("userFloat('dxyS')", float, doc="dxy/err (with sign) wrt first PV, in cm", precision=10),
        originalIndex = Var("userInt('originalIndex')", int, doc="index in the original Muon or Track collection", precision=10)
        ),
)


muonsTracksBParkSequence = cms.Sequence(muonsTracksBPark)
muonsTracksBParkTables = cms.Sequence(muonsTracksBParkTable)



