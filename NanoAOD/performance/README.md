## Performance testing
This folder contains a set of tools and instructions to measure 
peformance of NanoAOD production in various configurations

Use Bmm5/NanoAOD/performance/make_report.py to extract results.

## Results
Reference machine: vocms0118

### Monte Carlo

Dataset: BsToMuMu_BMuonFilter

#### Time per event for the event loop
| Tag           | Production Version | Reference NanoAOD | NanoAOD + Customizations |
| ------------- | ------------------ | ----------------- | ------------------------ |
| NanoAODv12-V01 |     526           |  0.096 sec/event  |       0.247 sec/event    |
| NanoAODv10-V03 |     523           |  0.053 sec/event  |       0.130 sec/event    |
| NanoAODv10-V02 |     522           |  0.039 sec/event  |  0.170 sec/event         |
| NanoAODv10-V01 |     521           |    |           |
| NanoAODv9-V03 |      519           |  0.136 sec/event  |       0.189 sec/event    |
| NanoAODv8-V02 |      516           |  0.078 sec/event  |       0.117 sec/event    |
| NanoAODv6-V18 |      515           |  0.073 sec/event  |       0.112 sec/event    |
| NanoAODv8-V01 |      514           |  0.074 sec/event  |       0.110 sec/event    |
| NanoAODv6-V17 |      513           |  0.070 sec/event  |       0.109 sec/event    |
| NanoAODv6-V14 |      511           |  0.070 sec/event  |       0.135 sec/event    |

#### File size per event 
| Tag           | Production Version | Reference NanoAOD | NanoAOD + Customizations |
| ------------- | ------------------ | ----------------- | ------------------------ | 
| NanoAODv12-V01 |     526           |   1.5 kB/event    |         3.6 kB/event     | 
| NanoAODv10-V03 |     523           |   1.3 kB/event    |         3.4 kB/event     | 
| NanoAODv10-V02 |     522           |   1.1 kB/event    |         3.2 kB/event     | 
| NanoAODv9-V03 |      519           |   1.2 kB/event    |         2.4 kB/event     | 
| NanoAODv6-V18 |      515           |                   |         1.9 kB/event     |
| NanoAODv8-V01 |      514           |                   |         2.0 kB/event     |
| NanoAODv6-V17 |      513           |                   |         1.9 kB/event     |
| NanoAODv6-V14 |      511           |   1.1 kB/event    |         1.7 kB/event     |

#### Memory Usage (RSS)
| Tag           | Production Version | Reference NanoAOD | NanoAOD + Customizations |
| ------------- | ------------------ | ----------------- | ------------------------ |
| NanoAODv12-V01 |     526           |      1359 kB      |         1589 kB          |
| NanoAODv10-V03 |     523           |      1728 kB      |         1960 kB          |
| NanoAODv10-V02 |     522           |      1179 kB      |         1546 kB          |
| NanoAODv9-V03 |      519           |      1763 kB      |         1945 kB          |
| NanoAODv8-V02 |      516           |                   |         1915 kB          |
| NanoAODv6-V18 |      515           |                   |         1893 kB          |
| NanoAODv8-V01 |      514           |                   |         1839 kB          |
| NanoAODv6-V17 |      513           |                   |         1886 kB          |
| NanoAODv6-V14 |      511           |      1630 kB      |         1893 kB          |

### Data: Charmonium Run2018D
#### Time per event for the event loop
| Tag           | Production Version | Reference NanoAOD | NanoAOD + Customizations |
| ------------- | ------------------ | ----------------- | ------------------------ |
| NanoAODv9-V03 |      519           |  0.067 sec/event  |       0.161 sec/event    |
| NanoAODv6-V18 |      515           |  0.076 sec/event  |       0.179 sec/event    |
| NanoAODv8-V01 |      514           |  0.100 sec/event  |       0.201 sec/event    |
| NanoAODv6-V17 |      513           |  0.076 sec/event  |       0.182 sec/event    |
| NanoAODv6-V14 |      511           |  0.073 sec/event  |       0.160 sec/event    |

#### File size per event
| Tag           | Production Version | Reference NanoAOD | NanoAOD + Customizations |
| ------------- | ------------------ | ----------------- | ------------------------ |
| NanoAODv9-V03 |      519           |   1.2 kB/event    |       3.7 kB/event       | 
| NanoAODv6-V18 |      515           |                   |       2.5 kB/event       |
| NanoAODv8-V01 |      514           |                   |       2.6 kB/event       |
| NanoAODv6-V17 |      513           |                   |       2.5 kB/event       |
| NanoAODv6-V14 |      511           |   1.0 kB/event    |       2.1 kB/event       |

#### Memory Usage (RSS)
| Tag           | Production Version | Reference NanoAOD | NanoAOD + Customizations |
| ------------- | ------------------ | ----------------- | ------------------------ |
| NanoAODv9-V03 |      519           |      1466 kB      |         1645 kB          |
| NanoAODv6-V18 |      515           |                   |         1785 kB          |
| NanoAODv8-V01 |      514           |                   |         1790 kB          |
| NanoAODv6-V17 |      513           |                   |         1794 kB          |
| NanoAODv6-V14 |      511           |      1665 kB      |         1799 kB          |


## Testing procedure
cmsDriver and processing commands
* Monte Carlo
   * BsToMuMu_BMuonFilter
      * Standard NanoAODv12 + all Bmm customization using Run3Summer22EEMiniAODv3 as input
         * ```cmsDriver.py RECO --conditions 130X_mcRun3_2022_realistic_v5 --customise Configuration/DataProcessing/Utils.addMonitoring --datatier NANOAOD --era Run3,run3_nanoAOD_124 --eventcontent NANOAODSIM --filein /store/user/dmytro/tmp/store+mc+Run3Summer22EEMiniAODv3+BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+MINIAODSIM+124X_mcRun3_2022_realistic_postEE_v1-v2+2820000+0096d5dd-88d3-46a0-a8cc-255a3090c71e.root --fileout file:BsToMuMu_BMuonFilter_NanoAODv12.root --nThreads 1 -n 10000 --no_exec --python_filename BsToMuMu_BMuonFilter_NanoAODv12.py --scenario pp --step NANO --mc --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeDileptonPlusX --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeV0ForMuonFake --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBmmMuonId --customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))" --customise Validation/Performance/TimeMemoryInfo.py --customise_commands="process.Timing.summaryOnly = cms.untracked.bool(True)"```
         * ```cmsRun BsToMuMu_BMuonFilter_NanoAODv12.py > & BsToMuMu_BMuonFilter_NanoAODv12.log```
         * ```python3 Bmm5/NanoAOD/performance/make_report.py BsToMuMu_BMuonFilter_NanoAODv12.log```
* Data:
   * Charmonium Run2018D
      * Standard NanoAODv9 + all Bmm customization
         * ```cmsDriver.py step1 --filein /store/user/dmytro/tmp/store+data+Run2018D+Charmonium+MINIAOD+UL2018_MiniAODv2-v1+240000+D9C795D0-EAC3-2A47-A631-E314B7AA9883.root --fileout file:Run2018D_NanoAOD_bmm.root --data --eventcontent NANOAOD --datatier NANOAOD --conditions 106X_dataRun2_v35 --step NANO --nThreads 1 --era Run2_2018,run2_nanoAOD_106Xv2 --python_filename Run2018D_NanoAOD_bmm.py --no_exec -n 10000 --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeDileptonPlusX --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeV0ForMuonFake --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBmmMuonId --customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))" --customise Validation/Performance/TimeMemoryInfo.py --customise_commands="process.Timing.summaryOnly = cms.untracked.bool(True)"```
         * ```cmsRun Run2018D_NanoAOD_bmm.py >& Run2018D_NanoAOD_bmm.log &```
         * ```python3 Bmm5/NanoAOD/performance/make_report.py Run2018D_NanoAOD_bmm.log```
