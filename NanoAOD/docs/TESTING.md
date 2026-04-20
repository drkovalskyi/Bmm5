# Bmm5/NanoAOD Testing

## Current State

The test infrastructure provides automated regression testing with
entry-by-entry content comparison, output validation, selection yield
checks, comparison plots, and performance extraction. There is no
CI/CD pipeline - tests are run manually.

### Existing Tests

**KinematicFit.cpp** (`test/KinematicFit.cpp`)

Standalone C++ program that tests kinematic vertex fitting for dimuon
+ photon events. Creates mock magnetic field, builds
KinematicParticles from concrete track parameters and covariance
matrices, runs vertex fit, and validates fit probability and chi2/NDF.
Built by scram as a binary. Run directly after `scram b`. No
assertions - returns 0 if fits succeed, requires manual inspection.

**runTrackerSeededPFCandAnalyzer_cfg.py** (`test/`)

CMSSW EDAnalyzer config. Processes MINIAODSIM/MINIAOD input and
analyzes TrackerSeededPFCandidates for muon track matching with
kinematic cuts (pT > 5, |eta| < 0.8). Run with `cmsRun`. Not a
test - a study tool.

**skim_fakes.py** (`test/`)

CMSSW config demonstrating event filtering using MuonFakeFilter. Run
with `cmsRun`. Supports local and DAS input. Not a test - a workflow
example.

### Existing Comparison Tools

**validation/rootdiff** - Python 2 / PyROOT CLI for branch-by-branch
ROOT file comparison. Features: auto-discover all TTrees, statistics
(min/max/mean/std), fast mode, branch pattern filter, tree selection,
branch rename mapping (JSON), tolerance. Limitation: Python 2 only.

**postprocess/compare_root_files.py** - Python 3 / uproot entry-by-
entry comparison. NaN-aware via numpy. Reports first mismatch per
branch. Limitation: hardcoded tree name, no CLI, no statistics, no
pattern filtering.

### Gaps

- No CI/CD pipeline
- No unit tests for utility classes (XGBooster, Displacement,
  CloseTrackInfo, Candidate)
- No test for Python configuration correctness
- KinematicFit.cpp has no pass/fail assertions
- Content comparison is only wired in at smoke level (single-threaded
  for deterministic event order). Nominal/release levels skip it
  because multi-thread cmsRun reorders events nondeterministically.
- `extract_performance.py` looks for a standalone `hostname` file as
  a legacy mechanism; hostname is already recorded in `test_info.json`
  so the companion file is no longer needed.

## Test Infrastructure

### Design Principles

A NanoAOD production run serves three purposes simultaneously:

1. **Smoke test** - the full chain runs without errors
2. **Content reference** (smoke level) - single-threaded production
   preserves event order, so the output ROOT file can be compared
   entry-by-entry against the reference to detect any content change
3. **Rough performance** - time/event and RSS from the log are
   sufficient to flag large regressions

Reference output and logs are stored on EOS. Performance comparisons
are only meaningful on the same host, so the hostname is recorded.

### Test Levels

Event counts are derived from per-sample nominal values (currently
500 for all samples) scaled by a level multiplier:

| Level       | Multiplier | Events (500 nominal) | Purpose                    |
|-------------|------------|----------------------|----------------------------|
| smoke       | 5%         | 25                   | Quick sanity check         |
| nominal     | 100%       | 500                  | Commit comparison          |
| pre-release | 20x        | 10,000               | Performance gate (1 sample)|
| release     | 20x        | 10,000               | Tagged release (all samples)|

### Scripts

All test scripts live in `test/`.

#### test_config.py

Central configuration for the test infrastructure. Defines:

- **Test levels** - multiplier table (smoke/nominal/release)
- **Paths** - EOS reference base, local output directory
- **Default sample** - `BuToJpsiK` (used for smoke/nominal when no
  sample is specified)
- **Customizations** - Bmm-specific and performance monitoring
  customizations applied to all samples
- **Sample definitions** - per-sample process parameters (type,
  era, conditions, input files, event count, optional selections
  and plot specs)
- **Git helpers** - functions to detect version tags and commit
  hashes from the Bmm5 repo

#### compare_root_files.py

Unified ROOT file comparison tool. Merges the capabilities of
`validation/rootdiff` (feature-rich CLI) and
`postprocess/compare_root_files.py` (Python 3 / uproot / exact
matching) into a single tool.

Usage:
```sh
python3 test/compare_root_files.py file1.root file2.root [options]
```

Features:
- Auto-discovers all TTrees in both files (Events, Runs,
  LuminosityBlocks, etc.)
- Entry-by-entry exact comparison (default mode) using numpy with
  NaN-aware matching
- Reports added/removed branches in both directions
- Reports first mismatch per differing branch with entry index and
  values

Options:
- `-v, --verbose` - statistics (min/max/mean/std) for differing
  branches
- `-f, --fast` - compare only entry counts and compressed sizes,
  skip data read
- `-p, --pattern REGEX` - analyze only branches matching pattern
- `-t, --tree NAME` - analyze specific tree only
- `-m, --map FILE` - JSON file mapping renamed branches
- `-i, --ignore` - suppress missing branch warnings
- `--tolerance FLOAT` - absolute tolerance for float comparison
  (default: exact)
- `--labels LABEL1 LABEL2` - custom labels for the two files in
  diff output

Output:
- Per-branch status: `[OK]`, `[DIFFERENT]`, `[ADDED]`, `[REMOVED]`
- Summary counts: matching, different, added, removed
- Exit code: 0 = identical, 1 = differences found

Implementation: Python 3, uproot + numpy. No ROOT dependency.

#### run_nanoaod.py

Runs a single NanoAOD production job for testing.

Usage:
```sh
python3 test/run_nanoaod.py -p SAMPLE [-l LEVEL] [-o OUTPUT] [-n NEVENTS]
python3 test/run_nanoaod.py --list
```

What it does:
1. Reads sample definition from `test_config.py`
2. Removes any stale outputs (`nanoaod_test.root`,
   `nanoaod_test.log`, `nanoaod_test_cfg.py`, `test_info.json`) up
   front so a failed run cannot be masked by a leftover ROOT file
   from a previous successful run.
3. Generates CMSSW config via `cmsDriver.py` with all Bmm
   customizations and performance monitoring
4. Runs `cmsRun`
5. Saves to output directory: `nanoaod_test.root`,
   `nanoaod_test.log`, `nanoaod_test_cfg.py`, `test_info.json`

Default level: nominal. Event count is computed from the sample's
base count scaled by the level multiplier, or overridden with `-n`.
`-n all` (or equivalently `-n -1`) processes every event in the
input file (cmsDriver EOF).
Smoke level uses `--nThreads 1` (single-threaded) to preserve event
order for deterministic entry-by-entry content comparison. Other
levels use `--nThreads 16`.

The `test_info.json` metadata file records: sample, type, era,
conditions, input files, event count, level, git version, git
commit, timestamp, and hostname.

#### check_output.py

Validates a NanoAOD output file for basic correctness.

Usage:
```sh
python3 test/check_output.py output.root
```

Checks:
- File exists, is non-empty, is valid ROOT
- Expected trees present (Events, Runs, LuminosityBlocks)
- Key Bmm branches exist (mm_kin_mass, mm_kin_vtx_prob,
  mm_kin_lxy, mm_kin_slxy, mm_kin_alpha, mm_kin_l3d, mm_kin_sl3d)
- Sanity: no NaN in mass fields, vtxProb in [0,1], lxy >= 0
- Reports event count

Exit code: 0 = all checks pass, 1 = failures found.

#### check_yields.py

Counts events passing selection cuts defined per sample.

Usage:
```sh
python3 test/check_yields.py output.root --selections 'JSON' [--ref reference.root]
```

Uses PyROOT `TTree::GetEntries(selection)` to count events. When a
reference file is provided, reports the difference. Selections are
defined in `test_config.py` per sample; all MC samples and Data have
selections configured.

#### make_plots.py

Produces overlay comparison histograms between current and reference
output.

Usage:
```sh
python3 test/make_plots.py current.root reference.root --plots 'JSON' [-o OUTPUT_DIR]
```

Uses PyROOT `TTree::Draw`. Reference histograms drawn in black,
current in magenta, line width 2. Stats boxes positioned
side-by-side. Plot specs are defined in `test_config.py` per sample;
most samples have plots configured.

#### extract_performance.py

Extracts performance metrics from a cmsRun log.

Usage:
```sh
python3 test/extract_performance.py logfile [--rootfile output.root] [--json] [--skip N]
```

Extracts:
- Hostname (from a companion `hostname` file or `--hostname` flag)
- Number of events processed
- Time per event (from "Begin processing" timestamps, excluding
  first N events for warmup; default N=1)
- Event loop time (CPU and Real per event from TimeReport Event
  Summary)
- Peak RSS
- File size per event (when `--rootfile` is given)
- Per-module timing from TimeReport Module Summary (all modules
  with >= 1 ms/event, sorted by time)
- DileptonPlusXProducer internal timing breakdown (ms/event per
  sub-function)

Output: human-readable by default, `--json` for structured output.

#### run_validation.py

Orchestrates the full validation workflow for one or more samples.

Usage:
```sh
python3 test/run_validation.py -p SAMPLE [-l LEVEL] [-r REF_TAG] [-o OUTPUT] [-n NEVENTS] [-s]
python3 test/run_validation.py --ref <baseline> --commit <target> -l smoke -p SAMPLE
```

The intended workflow is to run validation **before committing**:
modify code, `scram b`, then run validation. The reference is
built from HEAD (last commit) and current output is produced from
the compiled working tree, so the comparison shows exactly what
your uncommitted changes introduce.

**Two-commit mode**: Use `--commit` to build and test a specific
commit instead of the current working tree. Combined with `--ref`,
this compares two arbitrary commits without requiring either to be
checked out:

```sh
# Compare two commits (smoke level, entry-by-entry)
python3 test/run_validation.py --ref abc1234 --commit def5678 -l smoke -p BdToJpsiKShort

# Compare a commit against HEAD (nominal level)
python3 test/run_validation.py --commit def5678 -l nominal -p BdToJpsiKShort
```

Both sides are built in separate CMSSW areas (via
`build_reference_cmssw`). Output is cached under
`<base>/<level>/<commit_hash[:7]>/<sample>/`, so repeated runs
skip production if output already exists with matching event count.

Default behavior by level:
- **smoke/nominal**: runs `DEFAULT_SAMPLE` (`BuToJpsiK`),
  compares uncommitted working tree output against a reference
  auto-produced from HEAD by building a separate CMSSW area
- **pre-release**: runs `PRE_RELEASE_SAMPLE` (`BsToMuMu`) with
  10,000 events, compares against the latest tag. Reports detailed
  per-module timing, event loop time, and file size per event.
  Use this as a performance gate before a release.
- **release**: runs all samples with input files, compares against
  a tagged reference on EOS (latest Bmm5 tag or specified with `-r`)

Note: the reference CMSSW area uses the current working tree's test
scripts (`test_config.py`, `run_nanoaod.py`) with the old compiled
plugins. If test configuration and plugin code change together
(e.g., adding a new customization and the module it requires), the
reference build may fail. In that case, use release-level testing
against a pre-built EOS reference.

Steps:
1. **Produce current output** - always run first so we know how many
   events the current job processed. Uses `run_nanoaod.py`; skips
   re-production if the Bmm5 plugin library mtime matches the cached
   test_info.json. For `--commit` mode, also validates stored nevents
   matches the requested value; for the tag-reference path the check
   is skipped (see step 2b).
2. **Ensure reference** - if reference output exists on disk, use it
   tentatively; no length check. If not, build a CMSSW area at the
   reference commit and run production there.
2b. **Verify processed-event counts agree** - parse the `Events:`
   line (input events processed) from both current and reference
   logs via `extract_performance.py`. If they differ, the reference
   is force-rebuilt with `-n <current_events>` so the comparison is
   apples-to-apples. This defers the expensive ref rebuild until
   after we know current actually completed, avoiding wasted work.
3. **Validate output** - run `check_output.py` on current output
4. **Content comparison** (smoke only) - run `compare_root_files.py`
   entry-by-entry against the reference. Smoke uses single-threaded
   production to guarantee deterministic event order. Pass/fail.
5. **Selection yields** - run `check_yields.py` if sample has
   selections defined (compare against reference)
6. **Comparison plots** - run `make_plots.py` if sample has plot
   specs and reference exists
7. **Performance** - run `extract_performance.py` on current and
   reference logs

The `-n all` / `-n -1` override processes every input event; the
post-hoc check above handles the case where the cached reference was
built with a different `-n`.

Output: `validation_report.txt` with per-sample pass/fail results.
Exit code: 0 if all samples pass, 1 if any fail.

The `--save` flag copies current output to the EOS reference
directory for the current git version tag.

### Reference Management

References are stored on EOS (`/eos/cms/store/user/dmytro/tmp/bmm_test_references`).
Two comparison modes:

1. **Tag reference** (release level) - output from a fixed tag or
   the last release tag. Shows overall accumulated changes since the
   last release. Updated when `--save` is used.

2. **Commit reference** (smoke/nominal levels) - output from the
   last committed code (HEAD). Run validation before committing:
   the reference shows what HEAD produces, and the current output
   shows what your uncommitted changes produce. Auto-produced by
   `run_validation.py` via building a separate CMSSW environment
   at HEAD.

References are directories containing `nanoaod_test.root`,
`nanoaod_test.log`, and `test_info.json`.

### Test Data

Sample definitions live in `test_config.py`. Six samples are
configured:

| Sample         | Type | Era/Campaign                    | Events |
|----------------|------|---------------------------------|--------|
| BsToMuMu       | MC   | Run3Summer22EE, 124X            | 500    |
| BsToJPsiPhi    | MC   | RunIII2024Summer24, 140X        | 500    |
| BuToJpsiK      | MC   | RunIII2024Summer24, 140X        | 500    |
| BdToJpsiKstar  | MC   | RunIII2024Summer24, 140X        | 500    |
| BdToJpsiKShort | MC   | RunIII2024Summer24, 140X        | 500    |
| Data           | data | Run2024G ParkingDoubleMuonLowMass0 | 500 |

The default sample for smoke and nominal tests is **BuToJpsiK**,
which also has selection cuts and plot specs defined.

## Additional Tests (Lower Priority)

### Configuration Import Tests

Verify Python configs load without errors:

```sh
python3 -c "import Bmm5.NanoAOD.DileptonPlusX_cff"
python3 -c "import Bmm5.NanoAOD.BmmMuonId_cff"
python3 -c "import Bmm5.NanoAOD.BmmV0ForMuonFake_cff"
python3 -c "import Bmm5.NanoAOD.nano_cff"
```

Catches syntax errors and missing dependencies. Could be added to
`run_validation.py` as a pre-flight check.

### Unit Tests

Add scram-compatible C++ unit tests for core library classes. Each
test should be a standalone binary in `test/BuildFile.xml` returning
exit code 0 on success.

Priority targets:
1. **KinematicFitResult** - extend existing KinematicFit.cpp with
   assertions (mass in expected range, vtxProb in [0,1])
2. **XGBooster** - load model, set features, verify prediction
3. **Displacement** - verify compute_displacement() with known inputs

### Physics Validation (Manual)

The 86+ validation scripts in `validation/` cover physics-level
checks: efficiency, trigger turn-on, MVA ROC curves, background
composition, mass resolution. These require expert judgment and
should be run when physics-sensitive parameters change.

## Summary

| Component              | Script                   | Purpose                              |
|------------------------|--------------------------|--------------------------------------|
| Configuration          | test_config.py           | Sample definitions, levels, paths    |
| Production runner      | run_nanoaod.py           | Run NanoAOD job for one sample       |
| Output validator       | check_output.py          | Structure and sanity checks          |
| Content comparison     | compare_root_files.py    | Branch-by-branch ROOT diff           |
| Selection yields       | check_yields.py          | Event counts for selection cuts      |
| Comparison plots       | make_plots.py            | Overlay histograms (current vs ref)  |
| Performance extractor  | extract_performance.py   | Time/event, RSS from log             |
| Validation orchestrator| run_validation.py        | Orchestrate full test workflow        |

The primary workflow: after a code change, run `run_validation.py`
with the appropriate level. It will ensure a reference exists,
produce output, validate it, check yields, generate plots, and
report performance - giving a pass/fail answer with detailed
information about any changes.
