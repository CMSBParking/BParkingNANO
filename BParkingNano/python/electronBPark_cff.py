import FWCore.ParameterSet.Config as cms

lowptElectronsWithSeed = cms.EDProducer(
  'PATLowPtElectronSeedingEmbedder',
  src = cms.InputTag('FIXME'),
  biasedSeeding = cms.InputTag("lowPtGsfElectronSeedValueMaps","ptbiased","RECO"),
  unbiasedSeeding = cms.InputTag("lowPtGsfElectronSeedValueMaps","unbiased","RECO"),
)

lowptElectronsForAnalysis = cms.EDFilter(
  'PATElectronSelector',
  src = cms.InputTag("lowptElectronsWithSeed"),
  cut = cms.string(''),
  )

pfElectronsForAnalysis = cms.EDFilter(
  'PATElectronSelector',
  src = cms.InputTag("FIXME"),
  cut = cms.string(""),
  )

electronsForAnalysis = cms.EDProducer(
  'ElectronMerger',
  lowptSrc = cms.InputTag('lowptElectronsForAnalysis'),
  pfSrc    = cms.InputTag('pfElectronsForAnalysis'),
  drForCleaning = cms.double(0.00),
  dzForCleaning = cms.double(0.00),
  useGsfModeForP4 = cms.bool(True),
)

electrons = cms.Sequence(
  (
    lowptElectronsWithSeed *
    lowptElectronsForAnalysis +
    pfElectronsForAnalysis 
  ) *
  electronsForAnalysis
)
