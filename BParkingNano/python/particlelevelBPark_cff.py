import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.particlelevel_cff import *


particleLevelBParkSequence = cms.Sequence(mergedGenParticles + genParticles2HepMC + particleLevel)
particleLevelBParkTables = cms.Sequence(rivetLeptonTable)
