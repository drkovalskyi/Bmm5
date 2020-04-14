#!/bin/bash
echo "================= CMSRUN starting jobNum=$1 ====================" | tee -a job.log

source /cvmfs/cms.cern.ch/cmsset_default.sh

BASE=$PWD

cd CMSSW_10_2_18/src
eval `scram runtime -sh`

cd $BASE

n=4

sub_procs=`seq 0 $( expr $n - 1 )`
pids=()

for sub_proc in $sub_procs; do
    # make sub-process dir
    job_name="job-$sub_proc"
    if [ ! -d $job_name ] ; then 
	mkdir $job_name;
    fi
    lumi=$( expr $n '*' "$1" + "$sub_proc" )
    echo "================= CMSRUN starting sub-process $sub_proc ====================" | tee -a job.log
    cd $job_name && \
	cmsRun ../QCD_Pt-30toInf_BmmGenFilter_TuneCP5_13TeV_pythia8-RunIIFall18GS_cfg.py jobNum=$lumi > RunIIFall18GS.log 2>&1 && \
	cmsRun ../QCD_Pt-30toInf_BmmGenFilter_TuneCP5_13TeV_pythia8-RunIIAutumn18DRPremix-DIGI_cfg.py > RunIIAutumn18DRPremix-DIGI.log 2>&1 && \
	cmsRun ../QCD_Pt-30toInf_BmmGenFilter_TuneCP5_13TeV_pythia8-RunIIAutumn18DRPremix-RECO_cfg.py > RunIIAutumn18DRPremix-RECO.log 2>&1 && \
	cmsRun -e -j FrameworkJobReport.xml ../QCD_Pt-30toInf_BmmGenFilter_TuneCP5_13TeV_pythia8-RunIIAutumn18MiniAOD_cfg.py > RunIIAutumn18MiniAOD.log 2>&1 &
    pid=$!
    echo "Sub-process $sub_proc started with PID: $pid"
    pids+=($pid)
done

good_output=""

for sub_proc in $sub_procs; do
    pid=${pids[$sub_proc]}
    wait $pid
    if [ $? -eq 0 ]; then
        echo "SUCCESS - Job $item with $pid exited with a status of $?"
	if [ ! -z "$good_output" ]; then
	    good_output+=","
	fi
	good_output+="\"file:job-$sub_proc/QCD_Pt-30toInf_BmmGenFilter_TuneCP5_13TeV_pythia8-RunIIAutumn18MiniAOD.root\""
    else
        echo "FAILED - Job $item with $pid exited with a status of $?"
    fi
done

if [ ! -z "$good_output" ]; then

    echo "================= Merging output of successful jobs ====================" | tee -a job.log
    cat > merge.py <<EOF
import FWCore.ParameterSet.Config as cms

process = cms.Process("MERGE")

# Tell the process which files to use as the source
process.source = cms.Source ("PoolSource",
          fileNames = cms.untracked.vstring( $good_output )
)

process.maxEvents = cms.untracked.PSet(
            input = cms.untracked.int32 (-1)
)

process.Out = cms.OutputModule("PoolOutputModule",
         fileName = cms.untracked.string ("QCD_Pt-30toInf_BmmGenFilter_TuneCP5_13TeV_pythia8-RunIIAutumn18MiniAOD.root")
)

process.end = cms.EndPath(process.Out)
EOF
    cmsRun -e -j FrameworkJobReport.xml merge.py
else
    echo "FAILED to produce any useful output"
    sleep 100
    exit 1
fi

# echo "================= CMSRUN starting RunIIFall18GS ====================" | tee -a job.log
# cmsRun -j RunIIFall18GS.log QCD_Pt-30toInf_BmmGenFilter_TuneCP5_13TeV_pythia8-RunIIFall18GS_cfg.py jobNum=$1

# echo "================= CMSRUN starting RunIIAutumn18DRPremix-DIGI ====================" | tee -a job.log
# cmsRun -j RunIIAutumn18DRPremix-DIGI.log QCD_Pt-30toInf_BmmGenFilter_TuneCP5_13TeV_pythia8-RunIIAutumn18DRPremix-DIGI_cfg.py 

# echo "================= CMSRUN starting RunIIAutumn18DRPremix-RECO ====================" | tee -a job.log
# cmsRun -j RunIIAutumn18DRPremix-RECO.log QCD_Pt-30toInf_BmmGenFilter_TuneCP5_13TeV_pythia8-RunIIAutumn18DRPremix-RECO_cfg.py

# echo "================= CMSRUN starting RunIIAutumn18MiniAOD ====================" | tee -a job.log
# cmsRun QCD_Pt-30toInf_BmmGenFilter_TuneCP5_13TeV_pythia8-RunIIAutumn18MiniAOD_cfg.py 

# echo "================= CMSRUN finished ====================" | tee -a job.log
