from __future__ import print_function
import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.NanoAOD.globals_cff import *
from PhysicsTools.NanoAOD.nano_cff import *
from PhysicsTools.NanoAOD.vertices_cff import *
from PhysicsTools.NanoAOD.NanoAODEDMEventContent_cff import *

##for gen and trigger muon
from PhysicsTools.BParkingNano.genparticlesBPark_cff import *
from PhysicsTools.BParkingNano.particlelevelBPark_cff import *
from PhysicsTools.BParkingNano.triggerObjectsBPark_cff import *
from PhysicsTools.BParkingNano.RecoTrgMuonCand_cff import *

## filtered input collections
from PhysicsTools.BParkingNano.electronsBPark_cff import * 
from PhysicsTools.BParkingNano.muonsBPark_cff import * ##to be replaced with slimmed muons filtered wrt trigger muon




nanoSequenceOnlyFullSim = cms.Sequence(triggerObjectBParkTables + l1bits)

nanoSequence = cms.Sequence(nanoMetadata + 
                            vertexSequence +                             
                            globalTables + vertexTables + 
                            triggerObjectBParkTables + l1bits)

nanoSequenceMC = cms.Sequence(particleLevelBParkSequence + genParticleBParkSequence + 
                              muonBParkMC + electronBParkMC +
                              globalTablesMC + genWeightsTable + genParticleBParkTables + particleLevelBParkTables + lheInfoTable) 



def nanoAOD_customizeMC(process):
    process.nanoSequence = cms.Sequence( process.nanoSequence + nanoSequenceMC)
    return process

def nanoAOD_customizeMuonTriggerBPark(process):
    process.nanoSequence = cms.Sequence( process.nanoSequence + RecoTrgMuonCandSequence + RecoTrgMuonCandTables)
    return process

def nanoAOD_customizeElectronFilteredBPark(process):
    process.nanoSequence = cms.Sequence( process.nanoSequence + electronsBParkSequence + electronBParkTables)
    return process

