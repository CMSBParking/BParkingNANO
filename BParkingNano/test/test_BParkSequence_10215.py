import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
process = cms.Process('BParkNANO',eras.Run2_2018)


isMC = False
wantFullRECO = False

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('PhysicsTools.BParkingNano.nanoBPark_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)


globaltag = '102X_dataRun2_Sep2018Rereco_v1'
inputFiles = cms.untracked.vstring('/store/data/Run2018D/ParkingBPH5/MINIAOD/20Mar2019-v1/120000/0071842F-3A26-0D43-92F2-E2376273008E.root')
outputFileNANO = cms.untracked.string('testBParkNANO_data_10215.root')
outputFileFEVT = cms.untracked.string('testBParkFullEvt_data_10215.root')

if isMC:
    globaltag = '102X_upgrade2018_realistic_v15'
    inputFiles = cms.untracked.vstring('/store/user/bainbrid/lowpteleid/BuToKJpsi_Toee_MuFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/crab_lowpteleid/190328_152903/0000/step3_inMINIAODSIM_99.root')
    outputFileNANO = cms.untracked.string('testBParkNANO_mc_10215.root')
    outputFileFEVT = cms.untracked.string('testBParkFullEvt_mc_10215.root')


# Input source
process.source = cms.Source("PoolSource",
    fileNames = inputFiles,
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('test_data_10213 nevts:100'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition
process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-RECO'),
        filterName = cms.untracked.string('')
    ),
    fileName = outputFileFEVT,
    outputCommands = (cms.untracked.vstring('keep *')),
    splitLevel = cms.untracked.int32(0)
)

process.NANOAODoutput = cms.OutputModule("NanoAODOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(9),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('NANOAOD'),
        filterName = cms.untracked.string('')
    ),
    fileName = outputFileNANO,
    outputCommands = process.NANOAODEventContent.outputCommands
)


# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, globaltag, '')


# Path and EndPath definitions
process.nanoAOD_step = cms.Path(process.nanoSequence)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)
process.NANOAODoutput_step = cms.EndPath(process.NANOAODoutput)


# Schedule definition
process.schedule = cms.Schedule(process.nanoAOD_step,process.endjob_step, process.NANOAODoutput_step)
if wantFullRECO:
    process.schedule = cms.Schedule(process.nanoAOD_step,process.endjob_step, process.FEVTDEBUGHLToutput_step, process.NANOAODoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

from PhysicsTools.BParkingNano.nanoBPark_cff import *
process = nanoAOD_customizeMuonTriggerBPark(process)
process = nanoAOD_customizeElectronFilteredBPark(process)
process = nanoAOD_customizeTrackFilteredBPark(process)

# customisation of the process.
if isMC:
    from PhysicsTools.BParkingNano.nanoBPark_cff import nanoAOD_customizeMC
    process = nanoAOD_customizeMC(process)


process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
