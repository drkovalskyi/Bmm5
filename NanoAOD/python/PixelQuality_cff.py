import FWCore.ParameterSet.Config as cms

# NOTE: 
#    All instances of FlatTableProducers must end with Table in their
#    names so that their product match the keep patterns in the default
#    event content. Otherwise you need to modify outputCommands in
#    NanoAODEDMEventContent or provide a custom event content to the
#    output module

# readStuckTBM is set automatically from the GlobalTag by
# nanoAOD_customizeDileptonPlusX (True for Run-3 data, False otherwise) -- no
# manual override needed. The default here is False so the module never tries
# to read the stuckTBM tag unless something explicitly enables it: the tag is
# absent from every GlobalTag and its IOVs only start mid-2018 (run 315690), so
# reading it where it is not loaded/covered would abort the job. When False the
# masked columns get the -1 sentinel; the base/dead payload, the geometry-driven
# totals, and the per-muon HitPattern flags all still work on every era.
pixelQualityTable = cms.EDProducer(
    "PixelQualityTableProducer",
    readPixelQuality          = cms.bool(True),
    readStuckTBM              = cms.bool(False),
    pixelQualityStuckTBMLabel = cms.string("stuckTBM"),
    pixelQualityBaseLabel     = cms.string(""),
    verbosePixelQuality       = cms.untracked.bool(False),
)

pixelQualityMcTable = pixelQualityTable.clone( readPixelQuality = cms.bool(False) )

PixelQualityTables   = cms.Sequence(pixelQualityTable)
PixelQualityMcTables = cms.Sequence(pixelQualityMcTable)
