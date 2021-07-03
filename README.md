# Bmm5
## Bx to mumu analysis code based on MiniAOD input data in CMS experiment at CERN

Production configuration for Run NanoAODv8 with the new B-jet
reggression and Bmm5 NanoAODv8-V02 analysis code

## Build Instructions 
* scram p CMSSW CMSSW_10_6_19_patch2
* cd CMSSW_10_6_19_patch2/src/
* cmsenv
* git cms-addpkg PhysicsTools/NanoAOD
* git cms-merge-topic drkovalskyi:Bmm5-CMSSW_10_2_18-V01
* git clone git@github.com:drkovalskyi/Bmm5.git --branch NanoAODv8-V02
* scram setup Bmm5/NanoAOD/external-tools/rabit.xml
* scram setup Bmm5/NanoAOD/external-tools/xgboost.xml
* scram b -j 8

## cmsDriver Options
Here is a list required cmsDriver options to get analysis specific parts included
* --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBxToMuMu 
* --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeV0ForMuonFake
* --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBmmMuonId 

## Optional Filtering
If you want to add event filtering to the commands below you just need to modify the step option the following way
* Monte Carlo: --step NANO,FILTER:Bmm5/NanoAOD/BxToMuMuFilter_cff.BxToMuMuFilterSequenceMC
* Data: --step NANO,FILTER:Bmm5/NanoAOD/BxToMuMuFilter_cff.BxToMuMuFilterSequence

## Processing examples
* Configuration
  * Monte Carlo
    * **RunIISummer20UL18MiniAOD**
      * era: Run2_2018,run2_nanoAOD_106Xv1
      * conditions: 106X_upgrade2018_realistic_v15_L1v1
    * **RunIISummer20UL17MiniAOD**
      * era: Run2_2017,run2_nanoAOD_106Xv1
      * conditions: 106X_mc2017_realistic_v8
    * **RunIISummer20UL16MiniAOD**
      * era: Run2_2016,run2_nanoAOD_106Xv1
      * conditions: 106X_mcRun2_asymptotic_v15
    * **RunIISummer20UL16MiniAODAPV** - dynamic strip inefficiency
      * era: Run2_2016,run2_nanoAOD_106Xv1
      * conditions: 106X_mcRun2_asymptotic_preVFP_v9
  * Data
    * **Run2018**
      * era: Run2_2018,run2_nanoAOD_106Xv1
      * conditiongs: 106X_dataRun2_v32
    * **Run2017**
      * era: Run2_2017,run2_nanoAOD_106Xv1
      * conditiongs: 106X_dataRun2_v32
    * **Run2016**
      * era: Run2_2016,run2_nanoAOD_106Xv1
      * conditiongs: 106X_dataRun2_v32
* Multithreading
  * Production is using 4 cores
* Examples
  * Monte Carlo
    * cmsDriver.py step1 --filein /store/user/dmytro/tmp/store+mc+RunIISummer20UL18MiniAOD+BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+MINIAODSIM+106X_upgrade2018_realistic_v11_L1v1-v1+240000+A85F6114-1A37-2149-B06E-ABF1CEB9EC77.root --fileout file:RunIISummer20UL18MiniAOD.root --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --conditions 106X_upgrade2018_realistic_v15_L1v1 --step NANO --nThreads 1 --era Run2_2018,run2_nanoAOD_106Xv1 --python_filename RunIISummer20UL18MiniAOD.py --no_exec -n 10000 --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBxToMuMu --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeV0ForMuonFake --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBmmMuonId --customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))" --customise Validation/Performance/TimeMemoryInfo.py --customise_commands="process.Timing.summaryOnly = cms.untracked.bool(True)"
    * cmsDriver.py step1 --filein /store/user/dmytro/tmp/store+mc+RunIISummer20UL17MiniAOD+BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+MINIAODSIM+106X_mc2017_realistic_v6-v1+30000+597D662E-63D8-D846-BC74-B23813E75C58.root --fileout file:RunIISummer20UL17MiniAOD.root --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --conditions 106X_mc2017_realistic_v8 --step NANO --nThreads 1 --era Run2_2017,run2_nanoAOD_106Xv1 --python_filename RunIISummer20UL17MiniAOD.py --no_exec -n 10000 --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBxToMuMu --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeV0ForMuonFake --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBmmMuonId --customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))" --customise Validation/Performance/TimeMemoryInfo.py --customise_commands="process.Timing.summaryOnly = cms.untracked.bool(True)"
    * cmsDriver.py step1 --filein /store/user/dmytro/tmp/store+mc+RunIISummer20UL16MiniAOD+BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+MINIAODSIM+106X_mcRun2_asymptotic_v13-v1+100000+AD8066F7-1BFE-5A48-8AED-F558884E2E0B.root --fileout file:RunIISummer20UL16MiniAOD.root --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --conditions 106X_mcRun2_asymptotic_v15 --step NANO --nThreads 1 --era Run2_2016,run2_nanoAOD_106Xv1 --python_filename RunIISummer20UL16MiniAOD.py --no_exec -n 10000 --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBxToMuMu --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeV0ForMuonFake --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBmmMuonId --customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))" --customise Validation/Performance/TimeMemoryInfo.py --customise_commands="process.Timing.summaryOnly = cms.untracked.bool(True)"
    * cmsDriver.py step1 --filein /store/user/dmytro/tmp/store+mc+RunIISummer20UL16MiniAODAPV+BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+MINIAODSIM+106X_mcRun2_asymptotic_preVFP_v8-v1+240000+85E7B549-CF87-C74C-8810-10110216BFFD.root --fileout file:RunIISummer20UL16MiniAODAPV.root --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --conditions 106X_mcRun2_asymptotic_preVFP_v9 --step NANO --nThreads 1 --era Run2_2016,run2_nanoAOD_106Xv1 --python_filename RunIISummer20UL16MiniAODAPV.py --no_exec -n 10000 --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBxToMuMu --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeV0ForMuonFake --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBmmMuonId --customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))" --customise Validation/Performance/TimeMemoryInfo.py --customise_commands="process.Timing.summaryOnly = cms.untracked.bool(True)" 
  * Data
    * cmsDriver.py step1 --filein /store/user/dmytro/tmp/store+data+Run2018D+Charmonium+MINIAOD+12Nov2019_UL2018-v1+130000+4C06C314-A060-AB4D-82CD-1DAE9745CED0.root --fileout file:Run2018_NanoAOD.root --data --eventcontent NANOAOD --datatier NANOAOD --conditions 106X_dataRun2_v32 --step NANO --nThreads 1 --era Run2_2018,run2_nanoAOD_106Xv1 --python_filename Run2018_NanoAOD.py --no_exec -n 10000 --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBxToMuMu --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeV0ForMuonFake --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBmmMuonId --customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))" --customise Validation/Performance/TimeMemoryInfo.py --customise_commands="process.Timing.summaryOnly = cms.untracked.bool(True)"
    * cmsDriver.py step1 --filein /store/user/dmytro/tmp/store+data+Run2017D+Charmonium+MINIAOD+09Aug2019_UL2017-v1+30000+29D91DD3-1DBE-2940-8577-ADD115CC30E3.root --fileout file:Run2017B-F_NanoAOD.root --data --eventcontent NANOAOD --datatier NANOAOD --conditions 106X_dataRun2_v32 --step NANO --nThreads 1 --era Run2_2017,run2_nanoAOD_106Xv1 --python_filename Run2017B-F_NanoAOD.py --no_exec -n 10000 --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBxToMuMu --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeV0ForMuonFake --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBmmMuonId --customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))" --customise Validation/Performance/TimeMemoryInfo.py --customise_commands="process.Timing.summaryOnly = cms.untracked.bool(True)"
