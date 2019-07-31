import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *

def ufloat(expr, precision=-1):
    return Var('userFloat("%s")' % expr, float, precision = precision)

tracks = cms.EDFilter('TrackProducer',
       trgMuon    = cms.InputTag("slimmedMuons"),
       tracks     = cms.InputTag("packedPFCandidates"),
       lostTracks = cms.InputTag("lostTracks"),
       RunParameters = cms.PSet(
          #skips events that do not have certain obj
          SkipNoTrkEvt=cms.bool(False),
          TrkPtCut=cms.double(0.8),    TrkEtaCut=cms.double(2.5),
          TrkDzCut=cms.double(1.0),    TrkNormChiMin=cms.int32(-1),
          TrkNormChiMax=cms.int32(100),
                            ),
)


trackTable = cms.EDProducer(
    "SimpleCandidateFlatTableProducer",
    src = cms.InputTag("tracks:myTracks"),
    cut = cms.string(""), #we should not filter on cross linked collections
    name = cms.string("Track"),
    doc  = cms.string("Our packedCandidate table"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the muons
    variables = cms.PSet(
        CandVars,
         vx = Var("vx()", float, doc="X vertex",precision=10),
         vy = Var("vy()", float, doc="Y vertex",precision=10),
         vz = Var("vz()", float, doc="Z vertex",precision=10),
         Sdxy = Var("dxy()/dxyError()", float, doc="dxy significance",precision=10),
        ),
    externalVariables = cms.PSet(),
)

tables = cms.Sequence(
    tracks
    +trackTable 
)
