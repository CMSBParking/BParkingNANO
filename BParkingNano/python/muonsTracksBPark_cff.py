import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *

muonsTracksBPark = cms.EDProducer('MuonTrackMerger',
                                  muons  = cms.InputTag("muonTrgSelector:SelectedMuons"),
                                  muonTransientTracks = cms.InputTag("muonTrgSelector:SelectedTransientMuons"),
                                  tracks = cms.InputTag("tracksBPark:SelectedTracks"),
                                  trackTransientTracks = cms.InputTag("tracksBPark:SelectedTransientTracks"),
                                  muonSelection = cms.string(''),
                                  trackSelection = cms.string('(userInt("isPacked") == 1 && bestTrack().found() > 10 && bestTrack().found() < 25 && pt() < 4. && eta() < 1.5 && eta() > -1.5 )'),
                                  sortOutputCollections = cms.bool(False)
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
        isPacked = Var("userInt('isPacked')",int,doc="track from packedCandidate collection", precision=10),
        isPF = Var("userInt('isPF')",int,doc="from pfMuon collection", precision=10),
        dz = Var("userFloat('dz')",float,doc="dz (with sign) wrt first PV, in cm", precision=10),
        dxy = Var("userFloat('dxy')",float,doc="dxy (with sign) wrt first PV, in cm", precision=10),
        dzS = Var("userFloat('dzS')", float, doc="dz/err (with sign) wrt first PV, in cm", precision=10),
        dxyS = Var("userFloat('dxyS')", float, doc="dxy/err (with sign) wrt first PV, in cm", precision=10),
        DCASig=Var("userFloat('DCASig')", float,doc="significance of xy-distance of closest approach wrt beamspot", precision=10),
        #dEdXStrip=Var("userFloat('dEdXStrip')", float,doc="dE/dX from strips of associated isolated track"),
        #dEdXPixel=Var("userFloat('dEdXPixel')", float,doc="dE/dX from pixels of associated isolated track"),
        ),
)


muonsTracksBParkSequence = cms.Sequence(muonsTracksBPark)
muonsTracksBParkTables = cms.Sequence(muonsTracksBParkTable)



