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
from PhysicsTools.BParkingNano.muonsBPark_cff import * 

## filtered input collections
from PhysicsTools.BParkingNano.electronsBPark_cff import * 
from PhysicsTools.BParkingNano.tracksBPark_cff import *

## B collections
from PhysicsTools.BParkingNano.BToKLL_cff import *

nanoSequenceOnlyFullSim = cms.Sequence(triggerObjectBParkTables + l1bits)

nanoSequence = cms.Sequence(nanoMetadata + 
                            vertexSequence +           
                            globalTables + vertexTables + 
                            triggerObjectBParkTables + l1bits)

nanoSequenceMC = cms.Sequence(particleLevelBParkSequence + genParticleBParkSequence + 
                              globalTablesMC + genWeightsTable + genParticleBParkTables + particleLevelBParkTables + lheInfoTable) 



def nanoAOD_customizeMuonTriggerBPark(process):
    process.nanoSequence = cms.Sequence( process.nanoSequence + muonBParkSequence + muonTriggerMatchedTables + muonBParkTables)
    return process

def nanoAOD_customizeTrackFilteredBPark(process):
    process.nanoSequence = cms.Sequence( process.nanoSequence + tracksBParkSequence + tracksBParkTables)
    return process

def nanoAOD_customizeElectronFilteredBPark(process):
    process.nanoSequence = cms.Sequence( process.nanoSequence + electronsBParkSequence + electronBParkTables)
    return process

def nanoAOD_customizeBToKLL(process):
    process.nanoBKeeSequence   = cms.Sequence( BToKEESequence    + BToKeeTable   )
    process.nanoBKMuMuSequence = cms.Sequence( BToKMuMuSequence  + BToKmumuTable )
    return process

def nanoAOD_customizeMC(process):
    for ipath in process._Process__paths.items():
        path = process.paths[ipath[0]]
        try:
            index  = path.index(process.nanoBKeeSequence)
            path.insert(index+1, electronBParkMC)
        except:
            pass

        path.insert(0, nanoSequenceMC)
        path.insert(3, muonBParkMC+tracksBParkMC)

