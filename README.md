# Bmm5
## Bx to mumu analysis code based on MiniAOD input data in CMS experiment at CERN

Production configuration for Run NanoAODv6 with the new B-jet
reggression and Bmm5 NanoAODv8-V01 analysis code

## Build Instructions 
* scram p CMSSW CMSSW_10_6_20
* cd CMSSW_10_6_20/src/
* cmsenv
* git cms-addpkg PhysicsTools/NanoAOD
* git clone git@github.com:drkovalskyi/Bmm5.git --branch NanoAODv8-V01
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
      * conditions: 106X_mcRun2_asymptotic_v15
    * **RunIISummer20UL17MiniAOD**
      * era: Run2_2017,run2_nanoAOD_106Xv1
      * conditions: 106X_mcRun2_asymptotic_v15
    * **RunIISummer20UL16MiniAOD**
      * era: Run2_2016,run2_nanoAOD_106Xv1
      * conditions: 106X_mcRun2_asymptotic_v15
* Examples
  * Monte Carlo
    * cmsDriver.py step1 --filein /store/user/dmytro/tmp/store+mc+RunIISummer20UL18MiniAOD+BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+MINIAODSIM+106X_upgrade2018_realistic_v11_L1v1-v1+240000+A85F6114-1A37-2149-B06E-ABF1CEB9EC77.root --fileout file:RunIISummer20UL18MiniAOD.root --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --conditions 106X_mcRun2_asymptotic_v15 --step NANO --nThreads 1 --era Run2_2018,run2_nanoAOD_106Xv1 --python_filename RunIISummer20UL18MiniAOD.py --no_exec -n 10000 --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBxToMuMu --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeV0ForMuonFake --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBmmMuonId --customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))" --customise Validation/Performance/TimeMemoryInfo.py --customise_commands="process.Timing.summaryOnly = cms.untracked.bool(True)"
    * cmsDriver.py step1 --filein /store/user/dmytro/tmp/store+mc+RunIISummer20UL17MiniAOD+BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+MINIAODSIM+106X_mc2017_realistic_v6-v1+30000+597D662E-63D8-D846-BC74-B23813E75C58.root --fileout file:RunIISummer20UL17MiniAOD.root --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --conditions 106X_mcRun2_asymptotic_v15 --step NANO --nThreads 1 --era Run2_2017,run2_nanoAOD_106Xv1 --python_filename RunIISummer20UL17MiniAOD.py --no_exec -n 10000 --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBxToMuMu --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeV0ForMuonFake --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBmmMuonId --customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))" --customise Validation/Performance/TimeMemoryInfo.py --customise_commands="process.Timing.summaryOnly = cms.untracked.bool(True)"
    * cmsDriver.py step1 --filein /store/user/dmytro/tmp/store+mc+RunIISummer20UL16MiniAOD+BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+MINIAODSIM+106X_mcRun2_asymptotic_v13-v1+100000+AD8066F7-1BFE-5A48-8AED-F558884E2E0B.root --fileout file:RunIISummer20UL16MiniAOD.root --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --conditions 106X_mcRun2_asymptotic_v15 --step NANO --nThreads 1 --era Run2_2016,run2_nanoAOD_106Xv1 --python_filename RunIISummer20UL16MiniAOD.py --no_exec -n 10000 --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBxToMuMu --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeV0ForMuonFake --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBmmMuonId --customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))" --customise Validation/Performance/TimeMemoryInfo.py --customise_commands="process.Timing.summaryOnly = cms.untracked.bool(True)"
    * 
* Data (Nano1June2019)
  * **Run2016 (B-H)**
    * era: Run2_2016,run2_nanoAOD_94X2016
    * conditions: 102X_dataRun2_v12
    * cmsDriver.py step1 --filein root://cms-xrd-global.cern.ch//store/data/Run2016C/Charmonium/MINIAOD/17Jul2018-v1/20000/B802A313-FB8A-E811-96D5-0CC47AC52D78.root --fileout file:output.root --data --eventcontent NANOAOD --datatier NANOAOD --conditions 102X_dataRun2_v12 --step NANO --nThreads 2 --era Run2_2016,run2_nanoAOD_94X2016 --python_filename Run2016B-H_NanoAOD.py --no_exec -n 1000 --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBxToMuMu --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeV0ForMuonFake --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBmmMuonId --customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))"
  * **Run2017 (B-F)**
    * era: Run2_2017,run2_nanoAOD_94XMiniAODv2
    * conditions: 102X_dataRun2_v12
    * cmsDriver.py step1 --filein root://cms-xrd-global.cern.ch//store/data/Run2017D/Charmonium/MINIAOD/31Mar2018-v1/90000/9899D06E-B837-E811-8DBC-0025905C42A2.root --fileout file:output.root --data --eventcontent NANOAOD --datatier NANOAOD --conditions 102X_dataRun2_v12 --step NANO --nThreads 2 --era Run2_2017,run2_nanoAOD_94XMiniAODv2 --python_filename Run2017B-F_NanoAOD.py --no_exec -n 1000 --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBxToMuMu --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeV0ForMuonFake --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBmmMuonId --customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))"
  * **Run2018 (ABC)**
    * era: Run2_2018,run2_nanoAOD_102Xv1
    * conditions: 102X_dataRun2_v12
    * cmsDriver.py step1 --filein root://cms-xrd-global.cern.ch//store/data/Run2018C/Charmonium/MINIAOD/17Sep2018-v1/270000/8B6549C2-A288-0C41-9794-DF39E3F54CCD.root --fileout file:output.root --data --eventcontent NANOAOD --datatier NANOAOD --conditions 102X_dataRun2_v12 --step NANO --nThreads 2 --era Run2_2018,run2_nanoAOD_102Xv1 --python_filename Run2018ABC_NanoAOD.py --no_exec -n 1000 --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBxToMuMu --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeV0ForMuonFake --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBmmMuonId --customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))"
  * **Run2018 (prompt D)**
    * era: Run2_2018,run2_nanoAOD_102Xv1
    * conditions: 102X_dataRun2_Prompt_v15
    * cmsDriver.py step1 --filein root://cms-xrd-global.cern.ch//store/data/Run2018D/Charmonium/MINIAOD/PromptReco-v2/000/325/022/00000/0F526EF2-A897-C84D-9921-B8DFC60000EF.root --fileout file:output.root --data --eventcontent NANOAOD --datatier NANOAOD --conditions 102X_dataRun2_Prompt_v15 --step NANO --nThreads 2 --era Run2_2018,run2_nanoAOD_102Xv1 --python_filename Run2018D_NanoAOD.py --no_exec -n 1000 --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBxToMuMu --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeV0ForMuonFake --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBmmMuonId --customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))"
