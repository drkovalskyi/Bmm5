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

Monte Carlo processing based RunIIAutumn18NanoAODv5 campaign
* All events
  * cmsDriver.py step1 --filein /store/user/dmytro/tmp/store+mc+RunIIAutumn18MiniAOD+QCD_Pt-50to80_EMEnriched_TuneCP5_13TeV_pythia8+MINIAODSIM+102X_upgrade2018_realistic_v15-v1+120000+E5F6DFC8-65CF-2A41-B0BE-E82E041CB012.root --fileout file:output.root --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --conditions 102X_upgrade2018_realistic_v19 --step NANO --nThreads 2 --era Run2_2018,run2_nanoAOD_102Xv1 --python_filename nano_mc.py --no_exec -n 1000 --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBxToMuMu --customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))"
* Filtered 
  * cmsDriver.py step1 --filein /store/user/dmytro/tmp/store+mc+RunIIAutumn18MiniAOD+QCD_Pt-50to80_EMEnriched_TuneCP5_13TeV_pythia8+MINIAODSIM+102X_upgrade2018_realistic_v15-v1+120000+E5F6DFC8-65CF-2A41-B0BE-E82E041CB012.root --fileout file:output.root --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --conditions 102X_upgrade2018_realistic_v19 --step NANO,FILTER:Bmm5/NanoAOD/BxToMuMuFilter_cff.BxToMuMuFilterSequenceMC --nThreads 2 --era Run2_2018,run2_nanoAOD_102Xv1 --python_filename nano_mc.py --no_exec -n 1000 --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBxToMuMu --customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))"
