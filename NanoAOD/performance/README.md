## Performance testing
This folder contains a set of tools and instructions to measure 
peformance of NanoAOD production in various configurations

Use Bmm5/NanoAOD/performance/make_report.py to extract results.

## Results
### Monte Carlo: BsToMuMu_RunIIAutumn18NanoAODv6
#### Time per event for the event loop (sec)
| Tag           | Production Version | NanoAOD path | Bmm module | V0 module | MuonId module |
| ------------- | ------------------ | ------------ | ---------- | --------- | ------------- |
| NanoAODv6-V16 |      513           |     0.070    |    0.020   |    0.019  |    0.000      |
| NanoAODv6-V14 |      511           |     0.070    |    0.040   |    0.025  |    0.000      |

#### Size per event (KB)
| Tag           | Production Version | NanoAOD | NanoAOD + Bmm | NanoAOD + Bmm + V0 |
| ------------- | ------------------ | ------- |-------------- | ------------------ | 
| NanoAODv6-V16 |      513           |         |               |       1.9          |
| NanoAODv6-V14 |      511           |   1.1   |     1.6       |       1.7          |

#### Memory Usage (RSS KB)
| Tag           | Production Version | NanoAOD | NanoAOD + Bmm | NanoAOD + Bmm + V0 | NanoAOD + Bmm + V0 + Id |
| ------------- | ------------------ | ------- | ------------- | ------------------ | ----------------------- |
| NanoAODv6-V16 |      513           |         |               |                    |         1886            |
| NanoAODv6-V14 |      511           |   1630  |     1812      |      1903          |         1893            |

### Data: Charmonium Run2018D
#### Time per event for the event loop (sec)
| Tag           | Production Version | NanoAOD path | Bmm module | V0 module | MuonId module |
| ------------- | ------------------ | ------------ | ---------- | --------- | ------------- |
| NanoAODv6-V16 |      513           |     0.076    |    0.082   |    0.024  |    0.000      |
| NanoAODv6-V14 |      511           |     0.073    |    0.062   |    0.025  |    0.000      |

#### Size per event (KB)
| Tag           | Production Version | NanoAOD | NanoAOD + Bmm | NanoAOD + Bmm + V0 |
| ------------- | ------------------ | ------- |-------------- | ------------------ | 
| NanoAODv6-V16 |      513           |         |               |       2.5          |
| NanoAODv6-V14 |      511           |   1.0   |     2.0       |       2.1          |

#### Memory Usage (RSS KB)
| Tag           | Production Version | NanoAOD | NanoAOD + Bmm | NanoAOD + Bmm + V0 | NanoAOD + Bmm + V0 + Id |
| ------------- | ------------------ | ------- | ------------- | ------------------ | ----------------------- |
| NanoAODv6-V16 |      513           |         |               |                    |         1794            |
| NanoAODv6-V14 |      511           |   1665  |     1697      |      1761          |         1799            |


## Testing procedure
cmsDriver and processing commands
* Monte Carlo
   * BsToMuMu_RunIIAutumn18NanoAODv6
      * Standard NanoAOD + Bmm + V0 muon fakes + Muon Id
         * ```cmsDriver.py step1 --filein /store/group/phys_muon/dmytro/tmp/store+mc+RunIIAutumn18MiniAOD+BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+MINIAODSIM+102X_upgrade2018_realistic_v15-v1+80000+22D26CAC-CC86-E44F-A1C7-F2C5BA567CB1.root --fileout file:BsToMuMu_bmm_fakes_and_ids.root --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --conditions 102X_upgrade2018_realistic_v20 --step NANO --nThreads 1 --era Run2_2018,run2_nanoAOD_102Xv1 --python_filename BsToMuMu_RunIIAutumn18NanoAODv6_bmm_with_fakes_and_ids.py --no_exec -n 10000 --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBxToMuMu --customise Validation/Performance/TimeMemoryInfo.py --customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))" --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeV0ForMuonFake --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBmmMuonId --customise_commands="process.Timing.summaryOnly = cms.untracked.bool(True)"```
	 * ```cmsRun BsToMuMu_RunIIAutumn18NanoAODv6_bmm_with_fakes_and_ids.py >& BsToMuMu_RunIIAutumn18NanoAODv6_bmm_with_fakes_and_ids.log```
	 * ```python3 Bmm5/NanoAOD/performance/make_report.py BsToMuMu_RunIIAutumn18NanoAODv6_bmm_with_fakes_and_ids.log```
* Data:
   * Charmonium Run2018D
      * Standard NanoAOD + Bmm + V0 muon fakes
         * ```cmsDriver.py step1 --filein /store/user/dmytro/tmp/store+data+Run2018D+Charmonium+MINIAOD+PromptReco-v2+000+325+022+00000+0F526EF2-A897-C84D-9921-B8DFC60000EF.root --fileout file:Run2018D_NanoAOD_bmm_with_fakes.root --data --eventcontent NANOAOD --datatier NANOAOD --conditions 102X_dataRun2_Prompt_v15 --step NANO --nThreads 1 --era Run2_2018,run2_nanoAOD_102Xv1 --python_filename Run2018D_NanoAOD_bmm_with_fakes_and_ids.py --no_exec -n 10000 --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBxToMuMu --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeV0ForMuonFake --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBmmMuonId --customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))" --customise Validation/Performance/TimeMemoryInfo.py --customise_commands="process.Timing.summaryOnly = cms.untracked.bool(True)"```
         * ```cmsRun Run2018D_NanoAOD_bmm_with_fakes_and_ids.py >& Run2018D_NanoAOD_bmm_with_fakes_and_ids.log &```
	 * ```python3 Bmm5/NanoAOD/performance/make_report.py Run2018D_NanoAOD_bmm_with_fakes_and_ids.log```
