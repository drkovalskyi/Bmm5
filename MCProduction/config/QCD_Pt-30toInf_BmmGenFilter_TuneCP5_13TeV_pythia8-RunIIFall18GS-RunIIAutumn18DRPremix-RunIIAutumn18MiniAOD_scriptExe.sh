#!/bin/bash
echo "================= CMSRUN starting jobNum=$1 ====================" | tee -a job.log

source /cvmfs/cms.cern.ch/cmsset_default.sh

BASE=$PWD

cd CMSSW_10_2_18/src
eval `scram runtime -sh`

cd $BASE

echo "================= CMSRUN starting RunIIFall18GS ====================" | tee -a job.log
cmsRun -j RunIIFall18GS.log QCD_Pt-30toInf_BmmGenFilter_TuneCP5_13TeV_pythia8-RunIIFall18GS_cfg.py jobNum=$1

echo "================= CMSRUN starting RunIIAutumn18DRPremix-DIGI ====================" | tee -a job.log
cmsRun -j RunIIAutumn18DRPremix-DIGI.log QCD_Pt-30toInf_BmmGenFilter_TuneCP5_13TeV_pythia8-RunIIAutumn18DRPremix-DIGI_cfg.py 

echo "================= CMSRUN starting RunIIAutumn18DRPremix-RECO ====================" | tee -a job.log
cmsRun -j RunIIAutumn18DRPremix-RECO.log QCD_Pt-30toInf_BmmGenFilter_TuneCP5_13TeV_pythia8-RunIIAutumn18DRPremix-RECO_cfg.py

echo "================= CMSRUN starting RunIIAutumn18MiniAOD ====================" | tee -a job.log
cmsRun -e -j FrameworkJobReport.xml QCD_Pt-30toInf_BmmGenFilter_TuneCP5_13TeV_pythia8-RunIIAutumn18MiniAOD_cfg.py 

echo "================= CMSRUN finished ====================" | tee -a job.log
