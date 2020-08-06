## Performance testing
This folder contains a set of tools and instructions to measure 
peformance of NanoAOD production in various configurations

## Results
### Monte Carlo: BsToMuMu_RunIIAutumn18NanoAODv6
#### Average Time Per Event
| Tag           | Production Version | NanoAOD | NanoAOD + Bmm | NanoAOD + Bmm + V0 |
| ------------- | ------------------ | ------- |-------------- | ------------------ | 
| Master        |                    | 0.148   | 0.163         | 0.291              | 
| NanoAODv6-V07 | 507                | 0.150   | 0.161         | 0.270              |
| NanoAODv6-V05 | 505                | 0.148   | 0.159         | -                  |

#### Size per event (KiB)
| Tag           | Production Version | NanoAOD | NanoAOD + Bmm | NanoAOD + Bmm + V0 |
| ------------- | ------------------ | ------- |-------------- | ------------------ | 
| Master        | -                  | 3.48    | 4.02          | 4.03               | 
| NanoAODv6-V07 | 507                | 3.48    | 3.98          | 3.99               |
| NanoAODv6-V05 | 505                | 3.48    | 3.89          | -                  |

## Testing procedure
cmsDriver and processing commands
* Monte Carlo
   * BsToMuMu_RunIIAutumn18NanoAODv6
      * Standard NanoAOD
         * ```cmsDriver.py step1 --filein /store/group/phys_muon/dmytro/tmp/store+mc+RunIIAutumn18MiniAOD+BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+MINIAODSIM+102X_upgrade2018_realistic_v15-v1+80000+22D26CAC-CC86-E44F-A1C7-F2C5BA567CB1.root --fileout file:BsToMuMu.root --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --conditions 102X_upgrade2018_realistic_v20 --step NANO --nThreads 1 --era Run2_2018,run2_nanoAOD_102Xv1 --python_filename BsToMuMu_RunIIAutumn18NanoAODv6.py --no_exec -n 1000 --customise Validation/Performance/TimeMemoryInfo.py --customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))" --customise_commands="process.Timing.summaryOnly = cms.untracked.bool(True)"```
         * ```cmsRun BsToMuMu_RunIIAutumn18NanoAODv6.py > & BsToMuMu_RunIIAutumn18NanoAODv6.log &```
         * ```python3 make_report.py BsToMuMu_RunIIAutumn18NanoAODv6.log```
      * Standard NanoAOD + Bmm
         * ```cmsDriver.py step1 --filein /store/group/phys_muon/dmytro/tmp/store+mc+RunIIAutumn18MiniAOD+BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+MINIAODSIM+102X_upgrade2018_realistic_v15-v1+80000+22D26CAC-CC86-E44F-A1C7-F2C5BA567CB1.root --fileout file:BsToMuMu_RunIIAutumn18NanoAODv6_bmm.root --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --conditions 102X_upgrade2018_realistic_v20 --step NANO --nThreads 1 --era Run2_2018,run2_nanoAOD_102Xv1 --python_filename BsToMuMu_RunIIAutumn18NanoAODv6_bmm.py --no_exec -n 1000 --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBxToMuMu --customise Validation/Performance/TimeMemoryInfo.py --customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))" --customise_commands="process.Timing.summaryOnly = cms.untracked.bool(True)"```
         * ```cmsRun BsToMuMu_RunIIAutumn18NanoAODv6_bmm.py > & BsToMuMu_RunIIAutumn18NanoAODv6_bmm.log &```
         * ```python3 make_report.py BsToMuMu_RunIIAutumn18NanoAODv6_bmm.log```
      * Standard NanoAOD + Bmm + V0 muon fakes
         * ```cmsDriver.py step1 --filein /store/group/phys_muon/dmytro/tmp/store+mc+RunIIAutumn18MiniAOD+BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+MINIAODSIM+102X_upgrade2018_realistic_v15-v1+80000+22D26CAC-CC86-E44F-A1C7-F2C5BA567CB1.root --fileout file:BsToMuMu_bmm_fakes.root --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --conditions 102X_upgrade2018_realistic_v20 --step NANO --nThreads 1 --era Run2_2018,run2_nanoAOD_102Xv1 --python_filename BsToMuMu_RunIIAutumn18NanoAODv6_bmm_with_fakes.py --no_exec -n 1000 --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBxToMuMu --customise Validation/Performance/TimeMemoryInfo.py --customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))" --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeV0ForMuonFake --customise_commands="process.Timing.summaryOnly = cms.untracked.bool(True)"```
         * ```cmsRun BsToMuMu_RunIIAutumn18NanoAODv6_bmm_with_fakes.py > & BsToMuMu_RunIIAutumn18NanoAODv6_bmm_with_fakes.log &```
         * ```cmsRun BsToMuMu_RunIIAutumn18NanoAODv6_bmm_with_fakes.py > & BsToMuMu_RunIIAutumn18NanoAODv6_bmm_with_fakes.log &```

* Data:
   * Charmonium Run2018D
      * Standard NanoAOD
         *
      * Standard NanoAOD + Bmm
         * ```cmsDriver.py step1 --filein /store/data/Run2018D/Charmonium/MINIAOD/PromptReco-v2/000/325/022/00000/0F526EF2-A897-C84D-9921-B8DFC60000EF.root --fileout file:Run2018D_NanoAOD_bmm.root --data --eventcontent NANOAOD --datatier NANOAOD --conditions 102X_dataRun2_Prompt_v15 --step NANO --nThreads 1 --era Run2_2018,run2_nanoAOD_102Xv1 --python_filename Run2018D_NanoAOD_bmm.py --no_exec -n 1000 --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBxToMuMu --customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))" --customise Validation/Performance/TimeMemoryInfo.py --customise_commands="process.Timing.summaryOnly = cms.untracked.bool(True)"```
      * Standard NanoAOD + Bmm + V0 muon fakes
         * ```cmsDriver.py step1 --filein /store/data/Run2018D/Charmonium/MINIAOD/PromptReco-v2/000/325/022/00000/0F526EF2-A897-C84D-9921-B8DFC60000EF.root --fileout file:Run2018D_NanoAOD_bmm_with_fakes.root --data --eventcontent NANOAOD --datatier NANOAOD --conditions 102X_dataRun2_Prompt_v15 --step NANO --nThreads 1 --era Run2_2018,run2_nanoAOD_102Xv1 --python_filename Run2018D_NanoAOD_bmm_with_fakes.py --no_exec -n 1000 --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBxToMuMu --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeV0ForMuonFake --customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))" --customise Validation/Performance/TimeMemoryInfo.py --customise_commands="process.Timing.summaryOnly = cms.untracked.bool(True)"```
