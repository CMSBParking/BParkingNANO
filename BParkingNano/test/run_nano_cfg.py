from FWCore.ParameterSet.VarParsing import VarParsing
import FWCore.ParameterSet.Config as cms

options = VarParsing('python')

options.register('isMC', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run this on real data"
)
options.register('isMCunbiased', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "keep trigerless evts. Remove trigger bias"
)
options.register('globalTag', 'NOTSET',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Set global tag"
)
options.register('wantSummary', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run this on real data"
)
options.register('wantFullRECO', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run this on real data"
)
options.register('reportEvery', 500,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "report every N events"
)
options.register('skip', 0,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "skip first N events"
)

options.register('decay', [],
    VarParsing.multiplicity.list,
    VarParsing.varType.string,
    "specify decay channel(s) to run. options: kll, kstarll, none and only_trg"
)


options.setDefault('maxEvents', 100)
options.setDefault('tag', 'test')
options.parseArguments()


if options.isMCunbiased:
  options.isMC=True


globaltag = '102X_dataRun2_v11' if not options.isMC else '102X_upgrade2018_realistic_v15'
if options._beenSet['globalTag']:
    globaltag = options.globalTag
extension = {False : 'data', True : 'mc'}
outputFileNANO = cms.untracked.string('_'.join(['BParkNANO', extension[options.isMC], options.tag])+'.root')
outputFileFEVT = cms.untracked.string('_'.join(['BParkFullEvt', extension[options.isMC], options.tag])+'.root')
if not options.inputFiles:
    options.inputFiles = [
 #'/store/data/Run2018C/ParkingBPH3/MINIAOD/05May2019-v1/30000/655E367A-A896-D34A-8A23-24844C642BAF.root'
  #'/store/data/Run2018B/ParkingBPH6/MINIAOD/05May2019-v2/240000/89BF81DF-9C97-0144-A932-D67DDF827E20.root'
# '/store/data/Run2016E/Charmonium/MINIAOD/17Jul2018-v1/00000/6E943712-C68A-E811-813F-0026181D28F0.root'
#  '/store/data/Run2017F/Charmonium/MINIAOD/31Mar2018-v1/90000/EEBEF0B0-1737-E811-B8FF-0025904C66A4.root'
#  '/store/data/Run2017F/DoubleMuonLowMass/MINIAOD/31Mar2018-v1/00000/000B9885-0A38-E811-AE30-AC1F6B23C94A.root'
'/store/data/Run2016C/Charmonium/MINIAOD/17Jul2018-v1/40000/E8FEC4EC-D78B-E811-9F0E-00010100097E.root'
] if not options.isMC else \
                         [
#   "/store/mc/RunIIAutumn18MiniAOD/BuToKee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15_ext1-v2/270000/FE888274-B6BF-9A42-A81C-CD7C989FE423.root"
# "/store/mc/RunIIAutumn18MiniAOD/BuToKee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/60000/5EA6D799-A04A-1D42-A76B-6A564B28B7AF.root"
  "/store/mc/RunIIAutumn18MiniAOD/BuToKee_MufilterPt6_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v1/10000/CC7F1E87-6C4D-354F-9AA1-0E792BC9AF0B.root",
# "/store/mc/RunIIAutumn18MiniAOD/BuToK_ToMuMu_MuFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PUPoissonAve20_102X_upgrade2018_realistic_v15-v2/80000/FA2FDA05-CBB6-D84E-99A1-5E4A2FD6D696.root"
]
'''from os import listdir
for dr in listdir("/eos/cms/store/group/cmst3/user/gkaratha/GEN_PAT_BToKPsi2S_MuMu_HardQCDPtHut9.0_noMuFilter_v1.0"):
  options.inputFiles.append("/store/group/cmst3/user/gkaratha/GEN_PAT_BToKPsi2S_MuMu_HardQCDPtHut9.0_noMuFilter_v1.0/"+dr)'''

annotation = '%s nevts:%d' % (outputFileNANO, options.maxEvents)

from Configuration.StandardSequences.Eras import eras
process = cms.Process('BParkNANO',eras.Run2_2018)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('PhysicsTools.BParkingNano.nanoBPark_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

# Input source
process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles),
    secondaryFileNames = cms.untracked.vstring(),
    skipEvents=cms.untracked.uint32(options.skip),
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(options.wantSummary),
)

process.nanoMetadata.strings.tag = annotation
# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string(annotation),
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
    outputCommands = (cms.untracked.vstring('keep *',
                                            'drop *_*_SelectedTransient*_*',
                     )),
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
    outputCommands = cms.untracked.vstring(
     'drop *',
#      'keep *HLT*',
     "keep nanoaodFlatTable_*Table_*_*",     # event data
     "keep nanoaodUniqueString_nanoMetadata_*_*",   # basic metadata
    )

)


# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, globaltag, '')


from PhysicsTools.BParkingNano.nanoBPark_cff import *
#guard
if len(options.decay)==0:
   print "provide decay"
   exit()

if not options.isMCunbiased:
   process = nanoAOD_customizeMuonTriggerBPark(process)
else:
   process = nanoAOD_customizeMuonTriggerUnBiased(process)

if options.decay[0]!="only_trg":
   process = nanoAOD_customizeElectronFilteredBPark(process)
   process = nanoAOD_customizeTrackFilteredBPark(process)
process = nanoAOD_customizeTriggerBitsBPark(process)

if len(options.decay)==0 or "kll" in options.decay:
  process = nanoAOD_customizeBToKLL(process)
if len(options.decay)==0 or "kstarll" in options.decay:
  process = nanoAOD_customizeBToKstarEE(process)
  process = nanoAOD_customizeBToKstarMuMu(process)

process.schedule = cms.Schedule()
filters=[]

# Path and EndPath definitions
if "kll" in options.decay or len(options.decay)==0:
   if not options.isMCunbiased:
      process.nanoAOD_KMuMu_step = cms.Path(process.nanoSequence + process.nanoBKMuMuSequence + CountBToKmumu )
      process.nanoAOD_Kee_step   = cms.Path(process.nanoSequence + process.nanoBKeeSequence   + CountBToKee   )
   else:
      process.nanoAOD_KMuMu_step = cms.Path(process.nanoSequence + process.nanoBKMuMuSequence + CountBToKmumuUnBiased )
      process.nanoAOD_Kee_step   = cms.Path(process.nanoSequence + process.nanoBKeeSequence + CountBToKeeUnBiased )
#   process.nanoAOD_KMuMu_step = cms.Path(process.nanoSequence )
#   process.nanoAOD_Kee_step   = cms.Path(process.nanoSequence )
   process.schedule.insert(0, process.nanoAOD_KMuMu_step)
   process.schedule.insert(0, process.nanoAOD_Kee_step)
   filters.append('nanoAOD_KMuMu_step')
   filters.append('nanoAOD_Kee_step')

if "kstarll" in options.decay or len(options.decay)==0:
   process.nanoAOD_KstarMuMu_step = cms.Path(process.nanoSequence + process.KstarToKPiSequence + process.nanoBKstarMuMuSequence + CountBToKstarMuMu )
   process.nanoAOD_KstarEE_step  = cms.Path(process.nanoSequence+ process.KstarToKPiSequence + process.nanoBKstarEESequence + CountBToKstarEE  )
   process.schedule.insert(0, process.nanoAOD_KstarMuMu_step)
   process.schedule.insert(0, process.nanoAOD_KstarEE_step)
   filters.append('nanoAOD_KstarMuMu_step')
   filters.append('nanoAOD_KstarEE_step')

if options.decay[0]=="none" or options.decay[0]=="only_trg":
    process.nanoAOD_dryrun_step = cms.Path(process.nanoSequence)
    process.schedule.insert(0, process.nanoAOD_dryrun_step)
    filters.append('nanoAOD_dryrun_step')

# customisation of the process.
if options.isMC and not options.isMCunbiased:
    from PhysicsTools.BParkingNano.nanoBPark_cff import nanoAOD_customizeMC
    nanoAOD_customizeMC(process)

if options.isMCunbiased:
   process.gen_step = cms.Path( process.nanoSequence + process.nanoSequenceMC)
   process.schedule.insert(0, process.gen_step)
   filters.append('gen_step')
   

process.endjob_step = cms.EndPath(process.endOfProcess)
process.NANOAODoutput_step = cms.EndPath(process.NANOAODoutput)
process.schedule.insert(0,process.endjob_step)


if options.wantFullRECO:
    process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)
    process.schedule.insert(0, process.FEVTDEBUGHLToutput_step)

process.schedule.insert(0,process.NANOAODoutput_step)

from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

process.NANOAODoutput.SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring( filters )
)


### from https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/3287/1/1/1/1/1.html
process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))
process.NANOAODoutput.fakeNameForCrab=cms.untracked.bool(True)    

process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
