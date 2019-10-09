# Bmm5
## Bx to mumu analysis code based on MiniAOD input data in CMS experiment at CERN
Recommended release: CMSSW_10_2_15

Supported CMSSW releases: 10_2_X

Build Instructions 
* scram p CMSSW CMSSW_10_2_15
* cd CMSSW_10_2_15/src/
* cmsenv
* git clone git@github.com:drkovalskyi/Bmm5.git
* scram b -j 8

Reference campaign: Run2 NanoAODv5

Note: all examples below keep only useful events for Bmm5 analysis (loose di-muon filter). To disable filtering remove step FILTER from cmsDriver command.

* Monte Carlo
  * **RunIIAutumn18NanoAODv5**: 
    * cmsDriver.py step1 --filein /store/user/dmytro/tmp/store+mc+RunIIAutumn18MiniAOD+QCD_Pt-50to80_EMEnriched_TuneCP5_13TeV_pythia8+MINIAODSIM+102X_upgrade2018_realistic_v15-v1+120000+E5F6DFC8-65CF-2A41-B0BE-E82E041CB012.root --fileout file:output.root --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --conditions 102X_upgrade2018_realistic_v19 --step NANO,FILTER:Bmm5/NanoAOD/BxToMuMuFilter_cff.BxToMuMuFilterSequenceMC --nThreads 2 --era Run2_2018,run2_nanoAOD_102Xv1 --python_filename nano_RunIIAutumn18NanoAODv5.py --no_exec -n 1000 --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBxToMuMu --customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))"
* Data (Nano1June2019)
  * **Run2016 (B-H)**
    * era: Run2_2016,run2_nanoAOD_94X2016
    * conditions: 102X_dataRun2_v11
    * cmsDriver.py step1 --filein root://cms-xrd-global.cern.ch//store/data/Run2016C/Charmonium/MINIAOD/17Jul2018-v1/20000/B802A313-FB8A-E811-96D5-0CC47AC52D78.root --fileout file:output.root --data --eventcontent NANOAOD --datatier NANOAOD --conditions 102X_dataRun2_v11 --step NANO,FILTER:Bmm5/NanoAOD/BxToMuMuFilter_cff.BxToMuMuFilterSequence --nThreads 2 --era Run2_2016,run2_nanoAOD_94X2016 --python_filename nano_Run2016B-H.py --no_exec -n 1000 --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBxToMuMu --customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))"
  * **Run2017 (B-F)**
    * era: Run2_2017,run2_nanoAOD_94XMiniAODv2
    * conditions: 102X_dataRun2_v11
    * cmsDriver.py step1 --filein root://cms-xrd-global.cern.ch//store/data/Run2017D/Charmonium/MINIAOD/31Mar2018-v1/90000/9899D06E-B837-E811-8DBC-0025905C42A2.root --fileout file:output.root --data --eventcontent NANOAOD --datatier NANOAOD --conditions 102X_dataRun2_v11 --step NANO,FILTER:Bmm5/NanoAOD/BxToMuMuFilter_cff.BxToMuMuFilterSequence --nThreads 2 --era Run2_2017,run2_nanoAOD_94XMiniAODv2 --python_filename nano_Run2017B-F.py --no_exec -n 1000 --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBxToMuMu --customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))"
  * **Run2018 (ABC)**
    * era: Run2_2018,run2_nanoAOD_102Xv1
    * conditions: 102X_dataRun2_v11
    * cmsDriver.py step1 --filein root://cms-xrd-global.cern.ch//store/data/Run2018C/Charmonium/MINIAOD/17Sep2018-v1/270000/8B6549C2-A288-0C41-9794-DF39E3F54CCD.root --fileout file:output.root --data --eventcontent NANOAOD --datatier NANOAOD --conditions 102X_dataRun2_v11 --step NANO,FILTER:Bmm5/NanoAOD/BxToMuMuFilter_cff.BxToMuMuFilterSequence --nThreads 2 --era Run2_2018,run2_nanoAOD_102Xv1 --python_filename nano_Run2018ABC.py --no_exec -n 1000 --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBxToMuMu --customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))"
  * **Run2018 (prompt D)**
    * era: Run2_2018,run2_nanoAOD_102Xv1
    * conditions: 102X_dataRun2_Prompt_v14
    * cmsDriver.py step1 --filein root://cms-xrd-global.cern.ch//store/data/Run2018D/Charmonium/MINIAOD/PromptReco-v2/000/325/022/00000/0F526EF2-A897-C84D-9921-B8DFC60000EF.root --fileout file:output.root --data --eventcontent NANOAOD --datatier NANOAOD --conditions 102X_dataRun2_Prompt_v14 --step NANO,FILTER:Bmm5/NanoAOD/BxToMuMuFilter_cff.BxToMuMuFilterSequence --nThreads 2 --era Run2_2018,run2_nanoAOD_102Xv1 --python_filename nano_Run2018D.py --no_exec -n 1000 --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBxToMuMu --customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))"
