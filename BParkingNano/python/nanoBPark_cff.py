from __future__ import print_function
import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.BParkingNano.muonsBPark_cff import * 
from PhysicsTools.NanoAOD.globals_cff import *
from PhysicsTools.NanoAOD.nano_cff import *
from PhysicsTools.BParkingNano.genparticlesBPark_cff import *
from PhysicsTools.BParkingNano.particlelevelBPark_cff import *
from PhysicsTools.NanoAOD.vertices_cff import *
from PhysicsTools.BParkingNano.triggerObjectsBPark_cff import *
from PhysicsTools.BParkingNano.RecoTrgMuonCand_cff import *
from PhysicsTools.NanoAOD.NanoAODEDMEventContent_cff import *
##all following only cleaned wrt trigger muon
#need electrons pT > 1 and lowpt
#need tracks




nanoSequenceOnlyFullSim = cms.Sequence(triggerObjectBParkTables + l1bits)

nanoSequence = cms.Sequence(nanoMetadata + 
                            vertexSequence + 
                            muonBParkTables + 
                            globalTables + vertexTables + 
                            triggerObjectBParkTables + l1bits)

nanoSequenceMC = cms.Sequence(particleLevelBParkSequence + genParticleBParkSequence + 
                              muonBParkMC + 
                              globalTablesMC + genWeightsTable + genParticleBParkTables + particleLevelBParkTables + lheInfoTable) 



def nanoAOD_customizeMC(process):
    process.nanoSequence.insert(process.nanoSequence.index(triggerObjectBParkTables), cms.Sequence(process.nanoSequenceMC))
    return process

def nanoAOD_customizeMuonTriggerBPark(process):
    process.nanoSequence = cms.Sequence( process.nanoSequence + RecoTrgMuonCandSequence + RecoTrgMuonCandTables)
    return process

