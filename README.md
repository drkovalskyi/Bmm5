# Bmm5
## Bx to mumu analysis code based on MiniAOD input data in CMS experiment at CERN

Production configuration for Run NanoAODv6 with the new B-jet
reggression and Bmm5 NanoAODv6-V01 analysis code

## Build Instructions 
* scram p CMSSW CMSSW_10_2_18
* cd CMSSW_10_2_18/src/
* cmsenv
* git cms-init
* git remote add mithep git@github.com:MiT-HEP/cmssw.git
* git cms-addpkg PhysicsTools/NanoAOD
* git fetch mithep
* git checkout -t mithep/Bmm5-CMSSW_10_2_18
* git clone git@github.com:drkovalskyi/Bmm5.git --branch NanoAODv6-V01
* scram b -j 8

## Filtering
If you want to add event filtering to the commands below you just need to modify the step option the following way
* Monte Carlo: --step NANO,FILTER:Bmm5/NanoAOD/BxToMuMuFilter_cff.BxToMuMuFilterSequenceMC
* Data: --step NANO,FILTER:Bmm5/NanoAOD/BxToMuMuFilter_cff.BxToMuMuFilterSequence

## Processing examples
* Monte Carlo
  * **RunIIAutumn18NanoAODv5**: 
    * cmsDriver.py step1 --filein /store/user/dmytro/tmp/store+mc+RunIIAutumn18MiniAOD+QCD_Pt-50to80_EMEnriched_TuneCP5_13TeV_pythia8+MINIAODSIM+102X_upgrade2018_realistic_v15-v1+120000+E5F6DFC8-65CF-2A41-B0BE-E82E041CB012.root --fileout file:output.root --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --conditions 102X_upgrade2018_realistic_v19 --step NANO --nThreads 2 --era Run2_2018,run2_nanoAOD_102Xv1 --python_filename nano_RunIIAutumn18NanoAODv5.py --no_exec -n 1000 --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBxToMuMu --customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))"
* Data (Nano1June2019)
  * **Run2016 (B-H)**
    * era: Run2_2016,run2_nanoAOD_94X2016
    * conditions: 102X_dataRun2_v12
    * cmsDriver.py step1 --filein root://cms-xrd-global.cern.ch//store/data/Run2016C/Charmonium/MINIAOD/17Jul2018-v1/20000/B802A313-FB8A-E811-96D5-0CC47AC52D78.root --fileout file:output.root --data --eventcontent NANOAOD --datatier NANOAOD --conditions 102X_dataRun2_v12 --step NANO --nThreads 2 --era Run2_2016,run2_nanoAOD_94X2016 --python_filename Run2016B-H_NanoAOD.py --no_exec -n 1000 --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBxToMuMu --customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))"
  * **Run2017 (B-F)**
    * era: Run2_2017,run2_nanoAOD_94XMiniAODv2
    * conditions: 102X_dataRun2_v12
    * cmsDriver.py step1 --filein root://cms-xrd-global.cern.ch//store/data/Run2017D/Charmonium/MINIAOD/31Mar2018-v1/90000/9899D06E-B837-E811-8DBC-0025905C42A2.root --fileout file:output.root --data --eventcontent NANOAOD --datatier NANOAOD --conditions 102X_dataRun2_v12 --step NANO --nThreads 2 --era Run2_2017,run2_nanoAOD_94XMiniAODv2 --python_filename Run2017B-F_NanoAOD.py --no_exec -n 1000 --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBxToMuMu --customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))"
  * **Run2018 (ABC)**
    * era: Run2_2018,run2_nanoAOD_102Xv1
    * conditions: 102X_dataRun2_v12
    * cmsDriver.py step1 --filein root://cms-xrd-global.cern.ch//store/data/Run2018C/Charmonium/MINIAOD/17Sep2018-v1/270000/8B6549C2-A288-0C41-9794-DF39E3F54CCD.root --fileout file:output.root --data --eventcontent NANOAOD --datatier NANOAOD --conditions 102X_dataRun2_v12 --step NANO --nThreads 2 --era Run2_2018,run2_nanoAOD_102Xv1 --python_filename Run2018ABC_NanoAOD.py --no_exec -n 1000 --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBxToMuMu --customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))"
  * **Run2018 (prompt D)**
    * era: Run2_2018,run2_nanoAOD_102Xv1
    * conditions: 102X_dataRun2_Prompt_v15
    * cmsDriver.py step1 --filein root://cms-xrd-global.cern.ch//store/data/Run2018D/Charmonium/MINIAOD/PromptReco-v2/000/325/022/00000/0F526EF2-A897-C84D-9921-B8DFC60000EF.root --fileout file:output.root --data --eventcontent NANOAOD --datatier NANOAOD --conditions 102X_dataRun2_Prompt_v15 --step NANO --nThreads 2 --era Run2_2018,run2_nanoAOD_102Xv1 --python_filename Run2018D_NanoAOD.py --no_exec -n 1000 --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBxToMuMu --customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))"
