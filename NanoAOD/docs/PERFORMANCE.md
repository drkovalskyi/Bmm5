# Bmm5/NanoAOD Performance

## Evaluation Methods

Performance is measured using three metrics, extracted from CMSSW
cmsRun log output by `performance/make_report.py`.

### Metrics

**Time per event (event loop)** - Wall clock time for the event
processing loop divided by number of events. Extracted from the
`Total loop:` line in TimeReport output. The total time (including
initialization/finalization) is also reported; if the event loop
accounts for less than 90% of total time, a warning is printed
suggesting more events be processed for reliable results.

**File size per event** - Output ROOT file size divided by event
count. Measured in kB/event. Directly impacts storage requirements
for large-scale production (~1 TB total dataset).

**Memory usage (RSS)** - Peak resident set size from CMSSW
MemoryCheck. The script tracks the maximum RSS value across all
MemoryCheck log lines.

### Per-Module Breakdown

`make_report.py` also extracts per-module timing from the
`nanoAOD_step` path in the TimeReport. It reports:

- Total nanoAOD_step time per event
- Bmm-specific modules (matched by patterns: BxToMuMu,
  ForMuonFake, BmmMuonId, Dileptons) with percentage of path time
- Top 10 contributors by execution time

This breakdown identifies which modules dominate processing time.

## Benchmark Procedure

Reference machine: vocms0118 (CentOS7) / vocms118 (Alma9).
Single-threaded, 10000 events, IMT disabled. vocms118 was migrated
to new hardware on 2026-04-20; results produced on/after that date
are not directly comparable to earlier rows.

### Monte Carlo (BsToMuMu_BMuonFilter)

```sh
cmsDriver.py RECO \
  --conditions auto:phase1_2022_realistic_postEE \
  --customise Configuration/DataProcessing/Utils.addMonitoring \
  --datatier NANOAOD \
  --era Run3,run3_nanoAOD_124 \
  --eventcontent NANOAODSIM \
  --filein /store/user/dmytro/tmp/store+mc+Run3Summer22EEMiniAODv3+BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+MINIAODSIM+124X_mcRun3_2022_realistic_postEE_v1-v2+2820000+0096d5dd-88d3-46a0-a8cc-255a3090c71e.root \
  --fileout file:BsToMuMu_BMuonFilter_NanoAOD.root \
  --nThreads 1 -n 10000 --no_exec \
  --python_filename BsToMuMu_BMuonFilter_NanoAOD.py \
  --scenario pp --step NANO --mc \
  --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeDileptonPlusX \
  --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeV0ForMuonFake \
  --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBmmMuonId \
  --customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))" \
  --customise Validation/Performance/TimeMemoryInfo.py \
  --customise_commands="process.Timing.summaryOnly = cms.untracked.bool(True)" \
  --customise=Bmm5/NanoAOD/nano_cff.run3_nanoAOD_124

cmsRun BsToMuMu_BMuonFilter_NanoAOD.py >& BsToMuMu_BMuonFilter_NanoAOD.log
python3 Bmm5/NanoAOD/performance/make_report.py BsToMuMu_BMuonFilter_NanoAOD.log
```

### Data (Charmonium Run2018D)

```sh
cmsDriver.py step1 \
  --filein /store/user/dmytro/tmp/store+data+Run2018D+Charmonium+MINIAOD+UL2018_MiniAODv2-v1+240000+D9C795D0-EAC3-2A47-A631-E314B7AA9883.root \
  --fileout file:Run2018D_NanoAOD_bmm.root \
  --data --eventcontent NANOAOD --datatier NANOAOD \
  --conditions auto:run2_data \
  --step NANO --nThreads 1 \
  --era Run2_2018,run2_nanoAOD_106Xv2 \
  --python_filename Run2018D_NanoAOD_bmm.py \
  --no_exec -n 10000 \
  --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeDileptonPlusX \
  --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeV0ForMuonFake \
  --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBmmMuonId \
  --customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))" \
  --customise Validation/Performance/TimeMemoryInfo.py \
  --customise_commands="process.Timing.summaryOnly = cms.untracked.bool(True)"

cmsRun Run2018D_NanoAOD_bmm.py >& Run2018D_NanoAOD_bmm.log
python3 Bmm5/NanoAOD/performance/make_report.py Run2018D_NanoAOD_bmm.log
```

## Historical Results

Per-release numbers are maintained in
`Bmm5/NanoAOD/performance/README.md`. That file has two sections:

- **Results on vocms118 (new hardware, 2026-04-20+)** - NanoAODv15-V02
  vs. re-measured NanoAODv15-V01 reference with and without Bmm
  customizations.
- **Historical results (pre-migration vocms118 and earlier)** - the
  full per-tag tables going back to NanoAODv6-V14.

Refer to that file for the latest numbers. The trend analysis below
still applies to the pre-migration history.

## Trend Analysis

### Time per Event (MC, with customizations)

The processing time has grown significantly across versions:

- **v6 era** (511-515): ~0.11 sec/event - baseline
- **v9** (519): 0.189 sec/event - 1.7x baseline
- **v10** (522-523): 0.13-0.17 sec/event - some improvement
- **v12** (526-530): 0.25-0.32 sec/event - 2.3-2.9x baseline
- **v14** (531-534): 0.38-0.39 sec/event - 3.5x baseline
- **v15** (535): 0.57 sec/event - 5.2x baseline

The customization overhead (Bmm modules relative to reference
NanoAOD) has also grown. In v12-V01, customizations added 0.15
sec/event over a 0.096 base. By v15, reference NanoAOD timing is
not available, but the total is 0.57 sec/event.

### File Size (MC, with customizations)

- **v6 era**: 1.7-2.0 kB/event
- **v9**: 2.4 kB/event
- **v10**: 3.2-3.4 kB/event
- **v12**: 3.6-5.6 kB/event
- **v14**: 6.7-7.2 kB/event
- **v15**: 7.7 kB/event - 4.5x the v6 baseline

Growth comes from both upstream NanoAOD additions and Bmm-specific
output branches.

### Memory (MC, with customizations)

- **v6-v9 era**: ~1.9 GB
- **v12**: 1.6 GB (improvement)
- **v14**: 2.2 GB
- **v15**: 2.8 GB

### Data vs MC

Data processing (Charmonium Run2018D) shows a similar growth pattern
but with somewhat different absolute values. The v15 data processing
(0.503 sec/event, 11.5 kB/event) shows even more dramatic file size
growth than MC (11.5 vs 7.7 kB/event), likely due to higher
multiplicity in Charmonium data events.

## Observations for Optimization

1. **Time growth is partly upstream**: Reference NanoAOD itself slowed
   from 0.039 (v10-V02) to 0.096 (v12-V01) sec/event. The Bmm
   customization overhead grew on top of this.

2. **Per-module profiling is available**: `make_report.py` already
   identifies the top time consumers in the nanoAOD_step path. This
   should be run on the current version to identify specific modules
   to target.

3. **File size growth**: 4.5x growth from v6 to v15 suggests either
   new output branches being added or existing branches becoming
   larger. A branch-by-branch size comparison between versions would
   identify the sources.

4. **Memory growth**: 2.8 GB peak RSS in v15 is manageable but
   trending upward. Worth monitoring to prevent issues on
   resource-constrained batch nodes.

5. **Missing reference baselines**: Recent versions (v14, v15) lack
   reference NanoAOD measurements, making it impossible to separate
   upstream vs Bmm-specific contributions. Future benchmarks should
   always include a reference run.
