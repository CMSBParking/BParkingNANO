import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('python')

options.register('isData', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run this on real data"
)
options.register('wantSummary', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run this on real data"
)
options.register('reportEvery', 1000,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "report every N events"
)

options.setDefault('maxEvents', 100)
options.setDefault('inputFiles', ['/store/data/Run2018B/ParkingBPH4/MINIAOD/05May2019-v2/230000/6B5A24B1-0E6E-504B-8331-BD899EB60110.root'])
options.setDefault('outputFile', 'bparking_nano.root')
options.parseArguments()

globaltag = '102X_dataRun2_Sep2018Rereco_v1' if options.isData else'102X_upgrade2018_realistic_v15' 

process = cms.Process('ParkingNano')

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, globaltag, '')

process.maxEvents = cms.untracked.PSet( 
   input = cms.untracked.int32(options.maxEvents) 
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(options.wantSummary),
)

process.source = cms.Source(
   "PoolSource",
   # replace 'myfile.root' with the source file you want to use
   fileNames = cms.untracked.vstring(
      options.inputFiles
   ),
   secondaryFileNames=cms.untracked.vstring(
   ),
   inputCommands=cms.untracked.vstring(
      'keep *',
      'drop *_ctppsPixelClusters_*_*',      
   ),
)

#
# Here goes the ntuplising code with selection stuff
#
process.ntuplesSeq = cms.Sequence(
)

process.p = cms.Path(
   process.ntuplesSeq
   )

#
# NanoAOD output module
#
from PhysicsTools.NanoAOD.NanoAODEDMEventContent_cff import NANOAODEventContent
process.NANOAODoutput = cms.OutputModule(
   "NanoAODOutputModule",
   compressionAlgorithm = cms.untracked.string('LZMA'),
   compressionLevel = cms.untracked.int32(9),
   SelectEvents = cms.untracked.PSet(
      SelectEvents = cms.vstring('p') # Select only events passing path p
      ),
   dataset = cms.untracked.PSet(
      dataTier = cms.untracked.string('NANOAOD'),
        filterName = cms.untracked.string('')
   ),
   fileName = cms.untracked.string(options.outputFile),
   outputCommands = cms.untracked.vstring(
      'drop *',
      "keep nanoaodFlatTable_*Table_*_*",     # event data
      "keep nanoaodUniqueString_nanoMetadata_*_*",   # basic metadata
   ),
)

process.endjob = cms.EndPath(
   process.NANOAODoutput
)

process.schedule = cms.Schedule(process.p, process.endjob)   
