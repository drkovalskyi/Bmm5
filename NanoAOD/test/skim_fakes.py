import FWCore.ParameterSet.Config as cms
import re, subprocess
process = cms.Process('SKIM')
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("Bmm5.NanoAOD.MuonFakeFilter_cfi")
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

dataset = "/LambdaBToPPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM"

# Input source
process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring()
)

files = subprocess.check_output('dasgoclient -query "file dataset=%s"' % dataset, shell=True)
for f in files.split('\n'):
    if re.search('\S',f):
        process.source.fileNames.append(f)

output_name = re.sub('\/','+', re.sub('^\/', '', dataset))

process.outputPath = cms.EndPath()

process.output = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('file:%s.root' % output_name),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('filterPath'))                              
)

process.filterPath = cms.Path( process.muonFakeFilter )
process.output_step = cms.EndPath( process.output )
process.schedule = cms.Schedule( process.filterPath, process.output_step )

# Spit out filter efficiency at the end.                                                                                                                                         
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
