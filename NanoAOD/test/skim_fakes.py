import FWCore.ParameterSet.Config as cms
import re, subprocess
process = cms.Process('SKIM')
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("Bmm5.NanoAOD.MuonFakeFilter_cfi")
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

dataset = None
# dataset = "/LambdaBToPPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM"

# Input source
process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/user/dmytro/tmp/store+mc+RunIIAutumn18MiniAOD+QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_13TeV_pythia8+MINIAODSIM+102X_upgrade2018_realistic_v15-v3+60000+4D8AFF1E-98C3-464B-99F6-370070C6405D.root'
    )
)

output_name = "output.root"
if dataset:
    files = subprocess.check_output('dasgoclient -query "file dataset=%s"' % dataset, shell=True)
    for f in files.split('\n'):
        if re.search('\S',f):
            process.source.fileNames.append(f)

    output_name = re.sub('\/','+', re.sub('^\/', '', dataset)) + ".root"

process.outputPath = cms.EndPath()

process.output = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string(output_name),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('filterPath'))                              
)

process.filterPath = cms.Path( process.muonFakeFilter )
process.output_step = cms.EndPath( process.output )
process.schedule = cms.Schedule( process.filterPath, process.output_step )

# Spit out filter efficiency at the end.                                                                                                                                         
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
