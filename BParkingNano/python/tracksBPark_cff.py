import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *

tracksBPark = cms.EDProducer('TrackMerger',
                             beamSpot   = cms.InputTag("offlineBeamSpot"),
                             trgMuon    = cms.InputTag("muonTrgSelector:trgMuons"),
                             tracks     = cms.InputTag("packedPFCandidates"),
                             lostTracks = cms.InputTag("lostTracks"),
                             trkPtCut = cms.double(0.5),    
                             trkEtaCut = cms.double(2.5),
                             dzTrg_cleaning = cms.double(-1.), ##next step change to 1
                             drTrg_Cleaning = cms.double(-0.4),  ##next step move to 0.01
                             dcaSig = cms.double(-100000),
                             trkNormChiMin = cms.int32(-1),
                             trkNormChiMax = cms.int32(-1)
                            )


trackBParkTable = cms.EDProducer(
    "SimpleCompositeCandidateFlatTableProducer",
    src = cms.InputTag("tracksBPark:SelectedTracks"),
    cut = cms.string(""),
    name = cms.string("ProbeTracks"),
    doc  = cms.string("track collection probe side for BPark after basic selection"),
    singleton = cms.bool(False),
    extension = cms.bool(False), 
    variables = cms.PSet(
         CandVars,
        vx = Var("daughter(0).vx()", float, doc="x coordinate of vertex position, in cm", precision=10),
        vy = Var("daughter(0).vy()", float, doc="y coordinate of vertex position, in cm", precision=10),
        vz = Var("daughter(0).vz()", float, doc="z coordinate of vertex position, in cm", precision=10),
        isPacked = Var("userInt('isPacked')",int,doc="track from packedCandidate collection", precision=10),
        isLostTrk = Var("userInt('isLostTrk')",int,doc="track from lostTrack collection", precision=10),
        dz = Var("userFloat('dz')",float,doc="dz (with sign) wrt first PV, in cm", precision=10),
        dxy = Var("userFloat('dxy')",float,doc="dxy (with sign) wrt first PV, in cm", precision=10),
        dzS = Var("userFloat('dzS')", float, doc="dz/err (with sign) wrt first PV, in cm", precision=10),
        dxyS = Var("userFloat('dxyS')", float, doc="dxy/err (with sign) wrt first PV, in cm", precision=10),
        DCASig=Var("userFloat('DCASig')", float,doc="significance of xy-distance of closest approach wrt beamspot or trgMuon", precision=10),
        #dEdXStrip=Var("userFloat('dEdXStrip')", float,doc="dE/dX from strips of associated isolated track"),
        #dEdXPixel=Var("userFloat('dEdXPixel')", float,doc="dE/dX from pixels of associated isolated track"),
        ),
)


tracksBParkSequence = cms.Sequence(tracksBPark)
tracksBParkTables = cms.Sequence(trackBParkTable)



