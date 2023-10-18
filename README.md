# Bmm5
## B to mumu analysis code for data produced by the CMS experiment at CERN

The code is based on the **NanoAODv12** campaign for using Run3 Prompt Reco and MC MiniAOD. 
The setup produces standard NanoAOD output with additional branches
specific to the analysis. The samples can be used as a replacement for the standard
NanoAOD samples produced by the CMS Central Production infrastructure.

## Build Instructions 
* scram p CMSSW CMSSW_13_0_13
* cd CMSSW_13_0_13/src/
* cmsenv
* git clone git@github.com:drkovalskyi/Bmm5.git --branch NanoAODv12-V01
* scram b -j 8

## cmsDriver Options
Here is a list required cmsDriver options to get analysis specific parts included
* `--customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeDileptonPlusX`
* `--customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeV0ForMuonFake`
* `--customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBmmMuonId`

## Processing examples
* Configuration
  * Monte Carlo
    * **MiniAODv4**
      * era: Run3
      * conditions: 124X_mcRun3_2022_realistic_postEE_v1
    * **MiniAODv3**
      * era: Run3,run3_nanoAOD_124
      * conditions: 124X_mcRun3_2022_realistic_postEE_v1
      * Example:
      	* ```cmsDriver.py RECO --conditions 130X_mcRun3_2022_realistic_v5 --customise Configuration/DataProcessing/Utils.addMonitoring --datatier NANOAOD --era Run3,run3_nanoAOD_124 --eventcontent NANOAODSIM --filein /store/user/dmytro/tmp/store+mc+Run3Summer22EEMiniAODv3+BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+MINIAODSIM+124X_mcRun3_2022_realistic_postEE_v1-v2+2820000+0096d5dd-88d3-46a0-a8cc-255a3090c71e.root --fileout file:/tmp/dmytro/BsToMuMu_BMuonFilter_NanoAODv12_test.root --nThreads 4 -n "-1" --no_exec --python_filename BsToMuMu_BMuonFilter_NanoAODv12_test.py --scenario pp --step NANO --mc --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeDileptonPlusX --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeV0ForMuonFake --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBmmMuonId --customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))"```
  * Data
    * **Run2022**
      * era: Run3
      * conditiongs: 124X_dataRun3_Prompt_v4
* Multithreading
  * Production is using 4 cores, but anything up to 16 works well (not tested above that)
* Examples
  * Monte Carlo
    * ```cmsDriver.py RECO --conditions 124X_mcRun3_2022_realistic_postEE_v1 --customise Configuration/DataProcessing/Utils.addMonitoring --datatier NANOAOD --era Run3 --eventcontent NANOAODSIM --filein /store/user/dmytro/tmp/store+mc+Run3Summer22MiniAODv3+TTTo2J1L1Nu_CP5_13p6TeV_powheg-pythia8+MINIAODSIM+pilot_124X_mcRun3_2022_realistic_v11-v2+2530000+7a9f777f-2d9a-40b2-9d5d-1567ecc49975.root --fileout file:TTTo2J1L1Nu_NanoAODv10.root --nThreads 4 -n 1000 --no_exec --python_filename TTTo2J1L1Nu_NanoAODv10.py --scenario pp --step NANO --mc --customise PhysicsTools/NanoAOD/V10/nano_cff.nanoAOD_customizeV10 --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeDileptonPlusX --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeV0ForMuonFake --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBmmMuonId --customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))" --customise Validation/Performance/TimeMemoryInfo.py --customise_commands="process.Timing.summaryOnly = cms.untracked.bool(True)"```
  * Data
    * ```cmsDriver.py RECO --conditions 124X_dataRun3_Prompt_v4 --customise Configuration/DataProcessing/Utils.addMonitoring --datatier NANOAOD --era Run3 --eventcontent NANOAOD --filein /store/user/dmytro/tmp/store+data+Run2022F+ParkingDoubleMuonLowMass0+MINIAOD+PromptReco-v1+000+360+393+00000+014d212c-1429-4947-a250-8b471780ff2a.root --fileout file:Run2022F+ParkingDoubleMuonLowMass0-PromptNanoAODv10.root --nThreads 4 -n 1000 --no_exec --python_filename Run2022F+ParkingDoubleMuonLowMass0-PromptNanoAODv10.py --scenario pp --step NANO --data --customise PhysicsTools/NanoAOD/V10/nano_cff.nanoAOD_customizeV10 --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeDileptonPlusX --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeV0ForMuonFake --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBmmMuonId --customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))" --customise Validation/Performance/TimeMemoryInfo.py --customise_commands="process.Timing.summaryOnly = cms.untracked.bool(True)"```

## Optional Filtering
If you want to add event filtering to the commands below you just need to modify the step option the following way
* Monte Carlo: `--step NANO,FILTER:Bmm5/NanoAOD/BxToMuMuFilter_cff.BxToMuMuFilterSequenceMC`
* Data: `--step NANO,FILTER:Bmm5/NanoAOD/BxToMuMuFilter_cff.BxToMuMuFilterSequence`

