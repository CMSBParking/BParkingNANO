import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *

tracksBPark = cms.EDProducer('TrackMerger',
                             beamSpot   = cms.InputTag("offlineBeamSpot"),
                             trgMuon    = cms.InputTag("muonTrgSelector:trgMatched"),
                             tracks     = cms.InputTag("packedPFCandidates"),
                             lostTracks = cms.InputTag("lostTracks"),
                             trkPtCut = cms.double(0.5),    
                             trkEtaCut = cms.double(2.4),
                             dzTrg_cleaning = cms.double(1.),
                             drTrg_cleaning = cms.double(0.4),
                             dcaSig_probe = cms.double(1.5),
                             dcaSig_tag = cms.double(0.5),
                             trkNormChiMin = cms.int32(-1),
                             trkNormChiMax = cms.int32(100)
                            )


trackBParkTable = cms.EDProducer(
    "SimpleCompositeCandidateFlatTableProducer",
    src = cms.InputTag("tracksBPark:ProbeSide"),
    cut = cms.string(""), #we should not filter on cross linked collections
    name = cms.string("ProbeTracks"),
    doc  = cms.string("track collection probe side for BPark after basic selection"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the muons
    variables = cms.PSet(
        pt = Var("daughter(0).pt()",float,doc="pt", precision=10),
        eta = Var("daughter(0).eta()",float,doc="eta", precision=10),
        phi = Var("daughter(0).phi()",float,doc="phi", precision=10),
        charge = Var("daughter(0).charge()",float,doc="charge", precision=10),
        mass = Var("daughter(0).mass()",float,doc="mass", precision=10),
        pdgId = Var("daughter(0).pdgId()",int,doc="PF pdgID", precision=10),
        vx = Var("daughter(0).vx()", float, doc="x coordinate of vertex position, in cm", precision=10),
        vy = Var("daughter(0).vy()", float, doc="y coordinate of vertex position, in cm", precision=10),
        vz = Var("daughter(0).vz()", float, doc="z coordinate of vertex position, in cm", precision=10),
        dz = Var("userFloat('dz')",float,doc="dz (with sign) wrt first PV, in cm", precision=10),
        dxy = Var("userFloat('dxy')",float,doc="dxy (with sign) wrt first PV, in cm", precision=10),
        dzS = Var("userFloat('dzS')", float, doc="dz/err (with sign) wrt first PV, in cm", precision=10),
        dxyS = Var("userFloat('dxyS')", float, doc="dxy/err (with sign) wrt first PV, in cm", precision=10),
        DCASig=Var("userFloat('DCASig')", float,doc="significance of xy-distance of closest approach wrt beamspot or trgMuon"),
        #dEdXStrip=Var("userFloat('dEdXStrip')", float,doc="dE/dX from strips of associated isolated track"),
        #dEdXPixel=Var("userFloat('dEdXPixel')", float,doc="dE/dX from pixels of associated isolated track"),
        ),
)


trackTagSideBParkTable = trackBParkTable.clone(
    src = cms.InputTag("tracksBPark:TagSide"),
    name = cms.string("TagTracks"),
    doc  = cms.string("track collection tag side for BPark after basic selection"),
)

tracksBParkSequence = cms.Sequence(tracksBPark)
tracksBParkTables = cms.Sequence(trackBParkTable + trackTagSideBParkTable)
