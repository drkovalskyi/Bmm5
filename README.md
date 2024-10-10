# Bmm5
## B to mumu analysis code for data produced by the CMS experiment at CERN

The code is based on the **NanoAODv12** campaign for using Run3 Prompt Reco and MC MiniAOD. 
The setup produces standard NanoAOD output with additional branches
specific to the analysis. The samples can be used as a replacement for the standard
NanoAOD samples produced by the CMS Central Production infrastructure.

## Build Instructions 
* scram p CMSSW CMSSW_14_0_16
* cd CMSSW_14_0_16/src/
* cmsenv
* git clone git@github.com:drkovalskyi/Bmm5.git --branch NanoAODv12-V06
* scram b -j 8

## cmsDriver Options
Here is a list required cmsDriver options to get analysis specific parts included
* `--customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeDileptonPlusX`
* `--customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeV0ForMuonFake`
* `--customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBmmMuonId`

## Processing examples
* Configuration by MiniAOD campaigns
  * Monte Carlo
    * **Run3Summer22MiniAODv4** 
      * era: Run3
      * conditions: 130X_mcRun3_2022_realistic_v5
      * Example:
      	* ```cmsDriver.py RECO --conditions 130X_mcRun3_2022_realistic_v5 --datatier NANOAOD --era Run3 --eventcontent NANOAODSIM --filein /store/mc/Run3Summer22MiniAODv4/DstarToD0Pi_D0To2Mu_MuFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen/MINIAODSIM/130X_mcRun3_2022_realistic_v5-v1/30000/71ec4425-d76b-446d-9d89-a7b250c56568.root --fileout file:/tmp/dmytro/test.root --nThreads 4 -n 1000 --no_exec --python_filename test.py --scenario pp --step NANO --mc --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeDileptonPlusX --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeV0ForMuonFake --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBmmMuonId --customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))"```
    * **Run3Summer23BPixMiniAODv4** (Run3Summer23BPixNanoAODv12)
      * era: Run3_2023
      * conditions: 130X_mcRun3_2023_realistic_postBPix_v2
    * **Run3Summer23MiniAODv4** (Run3Summer23NanoAODv12)
      * era: Run3_2023
      * conditions: 130X_mcRun3_2023_realistic_v14
    * **MiniAODv3**
      * era: Run3,run3_nanoAOD_124
      * conditions: 130X_mcRun3_2022_realistic_v5
      * Example:
      	* ```cmsDriver.py RECO --conditions 130X_mcRun3_2022_realistic_v5 --datatier NANOAOD --era Run3,run3_nanoAOD_124 --eventcontent NANOAODSIM --filein /store/user/dmytro/tmp/store+mc+Run3Summer22EEMiniAODv3+BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+MINIAODSIM+124X_mcRun3_2022_realistic_postEE_v1-v2+2820000+0096d5dd-88d3-46a0-a8cc-255a3090c71e.root --fileout file:/tmp/dmytro/BsToMuMu_BMuonFilter_NanoAODv12_test.root --nThreads 4 -n "-1" --no_exec --python_filename BsToMuMu_BMuonFilter_NanoAODv12_test.py --scenario pp --step NANO --mc --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeDileptonPlusX --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeV0ForMuonFake --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBmmMuonId --customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))"```
    * **Run-2** - RunIISummer20UL MiniAODv2
      * era: Run2_2018,run2_nanoAOD_106Xv2
      * conditions: auto:phase1_2018_realistic
      * Example:
      	* ```cmsDriver.py RECO --conditions auto:phase1_2018_realistic --datatier NANOAOD --era Run2_2018,run2_nanoAOD_106Xv2 --eventcontent NANOAODSIM --filein /store/group/phys_muon/dmytro/tmp/store+mc+RunIISummer20UL18MiniAODv2+BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+MINIAODSIM+106X_upgrade2018_realistic_v16_L1v1-v2+120000+B73A72CF-5947-A14A-8686-92D34790C8F7.root --fileout file:BsToMuMu.root --nThreads 16 -n "-1" --no_exec --python_filename BsToMuMu.py --scenario pp --step NANO --mc --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeDileptonPlusX --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeV0ForMuonFake --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBmmMuonId --customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))"```
  * Data
    * **Run2023**
      * era: Run3
      * conditiongs: 130X_dataRun3_PromptAnalysis_v1
      * Example:
      	* ```cmsDriver.py RECO --conditions 130X_dataRun3_PromptAnalysis_v1 --datatier NANOAOD --era Run3 --eventcontent NANOAOD --filein /store/data/Run2023D/ParkingDoubleMuonLowMass0/MINIAOD/PromptReco-v1/000/370/580/00000/a0cf8056-01ce-48f5-9b7e-570e225b5aba.root --fileout file:/tmp/dmytro/test_data.root --nThreads 4 -n 10000 --no_exec --python_filename test_data.py --scenario pp --step NANO --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeDileptonPlusX --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeV0ForMuonFake --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBmmMuonId --customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))"```
    * **Run2022**
      * era: Run3,run3_nanoAOD_124
      * conditions: 130X_dataRun3_v2
      * Example:
      	* ```cmsDriver.py RECO --conditions 130X_dataRun3_v2 --datatier NANOAOD --era Run3,run3_nanoAOD_124 --eventcontent NANOAOD --filein file:/eos/cms/store/user/dmytro/tmp/store+data+Run2022C+ParkingDoubleMuonLowMass0+MINIAOD+PromptReco-v1+000+357+271+00000+ea64a9c2-6b1f-4744-b4ea-41aa0e3c3e1b.root --fileout file:/tmp/dmytro/test_data.root --nThreads 4 -n 10000 --no_exec --python_filename test_data.py --scenario pp --step NANO --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeDileptonPlusX --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeV0ForMuonFake --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBmmMuonId --customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))"```

## Optional Filtering
If you want to add event filtering to the commands below you just need to modify the step option the following way
* Monte Carlo: `--step NANO,FILTER:Bmm5/NanoAOD/BxToMuMuFilter_cff.BxToMuMuFilterSequenceMC`
* Data: `--step NANO,FILTER:Bmm5/NanoAOD/BxToMuMuFilter_cff.BxToMuMuFilterSequence`

