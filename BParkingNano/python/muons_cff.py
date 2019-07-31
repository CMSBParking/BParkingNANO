import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *

def ufloat(expr, precision=-1):
    return Var('userFloat("%s")' % expr, float, precision = precision)

muons = cms.EDFilter('MuonProducer',
       trgMuon    = cms.InputTag("slimmedMuons"),
       muons     = cms.InputTag("slimmedMuons"),
       vertices   = cms.InputTag("offlineSlimmedPrimaryVertices"),
       RunParameters = cms.PSet(
          #skips events that do not have certain obj
          SkipNoMuEvt=cms.bool(False),
          MuPtCut=cms.double(1.0),     MuEtaCut=cms.double(2.5),
          MuDzCut=cms.double(1.0),    MuSoftQ=cms.bool(False),
                            ),
)


muonTable = cms.EDProducer(
    "SimpleCandidateFlatTableProducer",
    src = cms.InputTag("muons:myMuons"),
    cut = cms.string(""), #we should not filter on cross linked collections
    name = cms.string("Muon"),
    doc  = cms.string("Our packedCandidate table"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the muons
    variables = cms.PSet(
        CandVars,
         vx = Var("vx()", float, doc="X vertex",precision=10),
         vy = Var("vy()", float, doc="Y vertex",precision=10),
         vz = Var("vz()", float, doc="Z vertex",precision=10),
       
        ),
    externalVariables = cms.PSet(),
)

tables = cms.Sequence(
    muons
    +muonTable 
)
