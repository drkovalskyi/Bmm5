## Performance testing
This folder contains a set of tools and instructions to measure 
peformance of NanoAOD production in various configurations

Use Bmm5/NanoAOD/performance/make_report.py to extract results.

## Results
### Monte Carlo
#### Time per event for the event loop
| Tag           | Production Version | Reference NanoAOD | NanoAOD + Customizations |
| ------------- | ------------------ | ----------------- | ------------------------ |
| NanoAODv9-V03 |                    |  0.136 sec/event  |       0.189 sec/event    |
| NanoAODv8-V02 |      516           |  0.078 sec/event  |       0.117 sec/event    |
| NanoAODv6-V18 |      515           |  0.073 sec/event  |       0.112 sec/event    |
| NanoAODv8-V01 |      514           |  0.074 sec/event  |       0.110 sec/event    |
| NanoAODv6-V17 |      513           |  0.070 sec/event  |       0.109 sec/event    |
| NanoAODv6-V14 |      511           |  0.070 sec/event  |       0.135 sec/event    |

#### File size per event 
| Tag           | Production Version | Reference NanoAOD | NanoAOD + Customizations |
| ------------- | ------------------ | ----------------- | ------------------------ | 
| NanoAODv9-V03 |                    |                   |         2.4 kB/event     | 
| NanoAODv6-V18 |      515           |                   |         1.9 kB/event     |
| NanoAODv8-V01 |      514           |                   |         2.0 kB/event     |
| NanoAODv6-V17 |      513           |                   |         1.9 kB/event     |
| NanoAODv6-V14 |      511           |   1.1 kB/event    |         1.7 kB/event     |

#### Memory Usage (RSS)
| Tag           | Production Version | Reference NanoAOD | NanoAOD + Customizations |
| ------------- | ------------------ | ----------------- | ------------------------ |
| NanoAODv9-V03 |                    |      1763 kB      |         1945 kB          |
| NanoAODv8-V02 |      516           |                   |         1915 kB          |
| NanoAODv6-V18 |      515           |                   |         1893 kB          |
| NanoAODv8-V01 |      514           |                   |         1839 kB          |
| NanoAODv6-V17 |      513           |                   |         1886 kB          |
| NanoAODv6-V14 |      511           |      1630 kB      |         1893 kB          |

### Data: Charmonium Run2018D
#### Time per event for the event loop (sec)
| Tag           | Production Version | NanoAOD path | Bmm module | V0 module | MuonId module |
| ------------- | ------------------ | ------------ | ---------- | --------- | ------------- |
| NanoAODv6-V18 |      515           |     0.076    |    0.080   |    0.023  |    0.000      |
| NanoAODv8-V01 |      514           |     0.100    |    0.079   |    0.022  |    0.000      |
| NanoAODv6-V17 |      513           |     0.076    |    0.082   |    0.024  |    0.000      |
| NanoAODv6-V14 |      511           |     0.073    |    0.062   |    0.025  |    0.000      |

#### Size per event (KB)
| Tag           | Production Version | NanoAOD | NanoAOD + Bmm | NanoAOD + Bmm + V0 |
| ------------- | ------------------ | ------- |-------------- | ------------------ | 
| NanoAODv6-V18 |      515           |         |               |       2.5          |
| NanoAODv8-V01 |      514           |         |               |       2.6          |
| NanoAODv6-V17 |      513           |         |               |       2.5          |
| NanoAODv6-V14 |      511           |   1.0   |     2.0       |       2.1          |

#### Memory Usage (RSS KB)
| Tag           | Production Version | NanoAOD | NanoAOD + Bmm | NanoAOD + Bmm + V0 | NanoAOD + Bmm + V0 + Id |
| ------------- | ------------------ | ------- | ------------- | ------------------ | ----------------------- |
| NanoAODv6-V18 |      515           |         |               |                    |         1785            |
| NanoAODv8-V01 |      514           |         |               |                    |         1790            |
| NanoAODv6-V17 |      513           |         |               |                    |         1794            |
| NanoAODv6-V14 |      511           |   1665  |     1697      |      1761          |         1799            |


## Testing procedure
cmsDriver and processing commands
* Monte Carlo
   * BsToMuMu_RunIIAutumn18NanoAODv9
      * Standard NanoAOD + Complete set of channels
         * ```cmsDriver.py step1 --filein /store/user/dmytro/tmp/store+mc+RunIISummer20UL18MiniAOD+BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+MINIAODSIM+106X_upgrade2018_realistic_v11_L1v1-v1+240000+A85F6114-1A37-2149-B06E-ABF1CEB9EC77.root --fileout file:BsToMuMu.root --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --conditions 106X_upgrade2018_realistic_v15_L1v1 --step NANO --nThreads 1 --era Run2_2018,run2_nanoAOD_106Xv1 --python_filename BsToMuMu_RunIIAutumn18NanoAODv9.py --no_exec -n 10000 --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeDileptonPlusX --customise Validation/Performance/TimeMemoryInfo.py --customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))" --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeV0ForMuonFake --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBmmMuonId --customise_commands="process.Timing.summaryOnly = cms.untracked.bool(True)"```
	 * ```cmsRun BsToMuMu_RunIIAutumn18NanoAODv9.py >& BsToMuMu_RunIIAutumn18NanoAODv9.log```
	 * ```python3 Bmm5/NanoAOD/performance/make_report.py BsToMuMu_RunIIAutumn18NanoAODv9.log```
* Data:
   * Charmonium Run2018D
      * Standard NanoAOD + Bmm + V0 muon fakes
         * ```cmsDriver.py step1 --filein /store/user/dmytro/tmp/store+data+Run2018D+Charmonium+MINIAOD+UL2018_MiniAODv2-v1+240000+D9C795D0-EAC3-2A47-A631-E314B7AA9883.root --fileout file:Run2018D_NanoAOD_bmm.root --data --eventcontent NANOAOD --datatier NANOAOD --conditions 106X_dataRun2_v35 --step NANO --nThreads 1 --era Run2_2018,run2_nanoAOD_106Xv2 --python_filename Run2018D_NanoAOD_bmm.py --no_exec -n 10000 --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeDileptonPlusX --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeV0ForMuonFake --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBmmMuonId --customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))" --customise Validation/Performance/TimeMemoryInfo.py --customise_commands="process.Timing.summaryOnly = cms.untracked.bool(True)"```
         * ```cmsRun Run2018D_NanoAOD_bmm.py >& Run2018D_NanoAOD_bmm.log &```
	 * ```python3 Bmm5/NanoAOD/performance/make_report.py Run2018D_NanoAOD_bmm.log```
