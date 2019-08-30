import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.genparticles_cff import *


# for BPHPark start with merged particles (pruned + packed),
# where pruned contain K* states, but not final states, 
# and packed contain final states (K pi).
# then you save also final states (granddaughters)
finalGenParticlesBPark = finalGenParticles.clone(
  src = cms.InputTag("mergedGenParticles"),
  select = cms.vstring(
	"drop *",
        "keep++ (abs(pdgId) == 511 || abs(pdgId) == 521)",  #keep all B0(=511) and B+/-(521) + their daughters and granddaughters
   )
)

genParticleBParkTable = genParticleTable.clone(
  src = cms.InputTag("finalGenParticlesBPark")
)

genParticleBParkSequence = cms.Sequence(finalGenParticlesBPark)
genParticleBParkTables = cms.Sequence(genParticleBParkTable)

