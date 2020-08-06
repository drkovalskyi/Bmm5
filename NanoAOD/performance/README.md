## Performance testing
This folder contains a set of tools and instructions to measure 
peformance of NanoAOD production in various configurations

## Test procedure
cmsDriver commands:
* Monte Carlo
   * BsToMuMu_RunIIAutumn18NanoAODv6
      * Standard NanoAOD
      * Standard NanoAOD + Bmm
         * ```cmsDriver.py step1 --filein /store/group/phys_muon/dmytro/tmp/store+mc+RunIIAutumn18MiniAOD+BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+MINIAODSIM+102X_upgrade2018_realistic_v15-v1+80000+22D26CAC-CC86-E44F-A1C7-F2C5BA567CB1.root --fileout file:BsToMuMu_RunIIAutumn18NanoAODv6_bmm.root --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --conditions 102X_upgrade2018_realistic_v20 --step NANO --nThreads 1 --era Run2_2018,run2_nanoAOD_102Xv1 --python_filename BsToMuMu_RunIIAutumn18NanoAODv6_bmm.py --no_exec -n 1000 --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBxToMuMu --customise Validation/Performance/TimeMemoryInfo.py --customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))" --customise_commands="process.Timing.summaryOnly = cms.untracked.bool(True)"```
* Data:

## Producing reports


## Results
