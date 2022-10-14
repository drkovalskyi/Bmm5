# Bmm5
## B to mumu analysis code for data produced by the CMS experiment at CERN

The code is based on the **NanoAODv9** campaign for using RunII UltraLegacy processing 
of **MiniAODv2**. The setup produces standard NanoAOD output with additional branches
specific to the analysis. The samples can be used as a replacement for the standard
NanoAOD samples produced by the CMS Central Production infrastructure.

## Build Instructions 
* scram p CMSSW CMSSW_10_6_26
* cd CMSSW_10_6_26/src/
* cmsenv
* git cms-addpkg PhysicsTools/NanoAOD
* git cms-merge-topic drkovalskyi:Bmm5-CMSSW_10_6_26
* git clone git@github.com:drkovalskyi/Bmm5.git --branch NanoAODv9-V03
* scram setup Bmm5/NanoAOD/external-tools/rabit.xml
* scram setup Bmm5/NanoAOD/external-tools/xgboost.xml
* scram b -j 8

## cmsDriver Options
Here is a list required cmsDriver options to get analysis specific parts included
* `--customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeDileptonPlusX`
* `--customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeV0ForMuonFake`
* `--customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBmmMuonId`

## Processing examples
* Configuration
  * Monte Carlo
    * **RunIISummer20UL18MiniAODv2**
      * era: Run2_2018,run2_nanoAOD_106Xv2
      * conditions: 106X_upgrade2018_realistic_v16_L1v1
    * **RunIISummer20UL17MiniAODv2**
      * era: Run2_2017,run2_nanoAOD_106Xv2
      * conditions: 106X_mc2017_realistic_v9
    * **RunIISummer20UL16MiniAODv2**
      * era: Run2_2016,run2_nanoAOD_106Xv2
      * conditions: 106X_mcRun2_asymptotic_v17
    * **RunIISummer20UL16MiniAODAPVv2** - dynamic strip inefficiency
      * era: Run2_2016_HIPM,run2_nanoAOD_106Xv2
      * conditions: 106X_mcRun2_asymptotic_preVFP_v11
  * Data
    * **Run2018**
      * era: Run2_2018,run2_nanoAOD_106Xv2
      * conditiongs: 106X_dataRun2_v35
    * **Run2017**
      * era: Run2_2017,run2_nanoAOD_106Xv2
      * conditiongs: 106X_dataRun2_v35
    * **Run2016**
      * era: Run2_2016,run2_nanoAOD_106Xv2
      * conditiongs: 106X_dataRun2_v35
* Multithreading
  * Production is using 4 cores
* Examples
  * Monte Carlo
    * ```cmsDriver.py step1 --filein /store/user/dmytro/tmp/store+mc+RunIISummer20UL18MiniAOD+BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+MINIAODSIM+106X_upgrade2018_realistic_v11_L1v1-v1+240000+A85F6114-1A37-2149-B06E-ABF1CEB9EC77.root --fileout file:BsToMuMu.root --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --conditions 106X_upgrade2018_realistic_v15_L1v1 --step NANO --nThreads 1 --era Run2_2018,run2_nanoAOD_106Xv1 --python_filename BsToMuMu_RunIIAutumn18NanoAODv9.py --no_exec -n 10000 --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeDileptonPlusX --customise Validation/Performance/TimeMemoryInfo.py --customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))" --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeV0ForMuonFake --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBmmMuonId --customise_commands="process.Timing.summaryOnly = cms.untracked.bool(True)"```
  * Data
    * ```cmsDriver.py step1 --filein /store/user/dmytro/tmp/store+data+Run2018D+Charmonium+MINIAOD+UL2018_MiniAODv2-v1+240000+D9C795D0-EAC3-2A47-A631-E314B7AA9883.root --fileout file:Run2018D_NanoAOD_bmm.root --data --eventcontent NANOAOD --datatier NANOAOD --conditions 106X_dataRun2_v35 --step NANO --nThreads 1 --era Run2_2018,run2_nanoAOD_106Xv2 --python_filename Run2018D_NanoAOD_bmm.py --no_exec -n 10000 --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeDileptonPlusX --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeV0ForMuonFake --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBmmMuonId --customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))" --customise Validation/Performance/TimeMemoryInfo.py --customise_commands="process.Timing.summaryOnly = cms.untracked.bool(True)"```

## Optional Filtering
If you want to add event filtering to the commands below you just need to modify the step option the following way
* Monte Carlo: `--step NANO,FILTER:Bmm5/NanoAOD/BxToMuMuFilter_cff.BxToMuMuFilterSequenceMC`
* Data: `--step NANO,FILTER:Bmm5/NanoAOD/BxToMuMuFilter_cff.BxToMuMuFilterSequence`

