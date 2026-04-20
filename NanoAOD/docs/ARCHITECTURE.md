# Bmm5/NanoAOD Architecture

## Overview

The Bmm5/NanoAOD package implements CMS flavor physics reconstruction
for B meson rare decay studies (Bs/Bd -> mu mu and related channels).
It extends standard CMS NanoAOD production with custom reconstruction
algorithms, kinematic fitting, MVA-based selection, and specialized
data formats for scouting streams.

The code is a CMSSW package built with `scram`. It consists of C++
EDM plugins for event reconstruction, a shared library of utility
classes, Python configuration files, and ML model data files.

## Directory Layout

```
Bmm5/NanoAOD/
  interface/        C++ headers (utility classes, data structures)
  src/              C++ implementations (shared library)
  plugins/          CMSSW EDM plugins (producers, filters, analyzers)
  python/           Python configuration for CMSSW modules
  data/             ML models (XGBoost, TMVA BDT, ONNX), certification JSONs
  test/             Standalone tests and helper configs
  performance/      Benchmark procedures and historical results
  validation/       Physics validation studies (efficiency, triggers, MVA)
  postprocess/      Flat ntuple production and distributed job system
  external-tools/   External library configs (XGBoost, RABIT)
  docs/             Documentation
```

## Build System

**Main library** (`BuildFile.xml`): Compiles `src/*.cc` and
`interface/*.h` into a shared library. Depends on xgboost, boost,
clhep, and CMSSW packages (DataFormats, TrackingTools, RecoVertex,
MagneticField).

**Plugins** (`plugins/BuildFile.xml`): Compiles EDM plugins. Additional
dependencies on PhysicsTools/TensorFlow, DataFormats/NanoAOD,
DataFormats/Scouting, PhysicsTools/PatAlgos.

**Tests** (`test/BuildFile.xml`): Builds the `KinematicFit` test binary.

## Core Library (interface/ + src/)

### Candidate (Candidate.h/cc)

`bmm::Candidate` - unified particle container inheriting from
`reco::LeafCandidate`. Wraps PAT muons, scouting muons, PAT
electrons, packed PF candidates (hadrons), and bare reco tracks
through a common interface. Provides `index()`, `track()`,
`bestTrack()`, `genParticle()` accessors. Used throughout the
plugins to handle different particle types uniformly.

### Kinematic Fitting (KinFitUtils.h/cc, KinematicFitResult.h/cc)

`KinFitUtils.h` defines particle mass constants, matrix typedefs,
`KalmanVertexFitResult`, `DisplacementInformationIn3D`, and helper
functions (`jacobianSphToCart`, `getFTS`, `build_particle`,
`makeLorentzVectorFromPxPyPzM`).

`KinematicFitResult` stores fit output: validity, mass, refitted
daughter momenta, chi2/ndof, vertex position/error, Lxy with error
and significance, alpha angle relative to beam spot. The
`postprocess()` method computes displacement properties from the
fitted vertex.

`set_tree()` additionally validates the fit via `hasFiniteValues()`:
vertex position/error, mother momentum/mass, the full 7x7 kinematic
covariance, and every daughter's momentum and mass are scanned for
NaN/Inf. Fits that "converge" with non-finite state (seen with
ill-conditioned matrices) are rejected - the tree is cleared and the
result is marked invalid. `valid()` alone is not sufficient; use the
full path `set_tree()` -> `valid()`.

### Displacement (Displacement.h/cc)

`bmm::Displacement` wraps displacement information (decay length,
distance of closest approach, longitudinal IP, decay time, alpha
angle) with errors and significances. `bmm::Displacements` is a
named collection with lookup by vertex type.

`compute_displacement()` is the core algorithm: computes all
displacement quantities from secondary vertex to production vertex
with full error propagation.

### Close Track Analysis (CommonTools.h/cc)

`CloseTrack` stores per-track information near a vertex (DOCA to SV
and PV, vertex probability, beam spot IP significance).
`CloseTrackInfo` aggregates tracks and provides counting methods:
`nTracksByVertexProbability()`, `nTracksByDisplacementSignificance()`,
`nTracksByBetterMatch()`, `minDoca()`. These feed isolation variables
into the candidate.

Additional utilities: `dr_match()` for delta-R matching,
`get_mother()` / `find_common_ancestor()` for gen particle
navigation, `get_pixel_pattern()` for hit pattern extraction,
`fill_track_info()` for track property storage.

### XGBoost Interface (XGBooster.h/cc)

Wraps the XGBoost C API for inference. Loads a model file and
optional feature name file (JSON). Features are set by name via
`set()`, `predict()` runs inference and auto-resets values to NaN
for safety. Used by BmmMuonIdProducer and DileptonPlusXProducer.

### Scouting Data Conversion (ScoutingDataHandling.h/cc)

Converts between scouting and standard CMSSW formats:
`makeRecoTrack(Run3ScoutingTrack/Muon)`, `makePatMuon()`,
`makeRecoVertex()`. Enables the scouting pipeline to reuse the same
reconstruction algorithms as the offline pipeline.

### Track Occupancy (TrackOccupancyInfo.h)

`TrackOccupancyData` stores per-event track counts by quality
criteria. `TrackOccupancyMetaData` defines the eta-phi binning grid.
Used by TrackOccupancyLumiProducer for monitoring.

## Plugins

### Main Reconstruction

**DileptonPlusXProducer** (EDProducer) - Primary reconstruction
plugin. Reconstructs Bs/Bd -> mu mu + X decay chains:

1. Selects muons/electrons passing quality cuts
2. Forms dilepton pairs and performs kinematic vertex fit
3. Optionally adds hadrons (K, pi) or photons for extended decay modes
4. Computes displacement metrics (Lxy, alpha, decay time)
5. Calculates isolation from close track analysis
6. Runs BDT (TMVA) and XGBoost MVA for signal/background separation
7. Matches to gen particles (MC mode)
8. Outputs `pat::CompositeCandidate` collections and NanoAOD FlatTables

Consumes: slimmedMuons, slimmedElectrons, slimmedPhotons,
packedPFCandidates, primaryVertices, beamSpot, genParticles, trigger
results, L1 muons.

Produces collections for: DiMuons, DiMuonMc, DiMuonPhi, DiMuonKs,
DiMuonPhoton, BToJpsiKs, and others.

Timing instrumentation is gated by the `BMM_PROFILING` compile flag.
Without the flag (default production build) all `perfTimer(...)` calls
and `PerfClock::now()` calls collapse to no-ops under the optimizer.
With `-DBMM_PROFILING` the plugin emits a single-line per-stream
timing report to stderr at `endStream()` listing per-phase totals
(`fit_kin`, `fit_jpsikk`, `fit_phill`, `fit_jpsiks`, `fillBtoKllInfo`,
`fillBtoLLhhInfo`, `buildLLXCandidates`, `mmGamma`, `dstar`, `kstar`,
`mmm`, `isolation`, `tnp`, `ee`, `emu`, `hh`). Build with
`scram b -j 16 USER_CXXFLAGS=-DBMM_PROFILING`.

**ScoutingDileptonPlusXProducer** (EDProducer) - Parallel
implementation for scouting data. Reads Run3Scouting collections,
converts via ScoutingDataHandling, then applies the same
reconstruction logic as DileptonPlusXProducer.

**BmmV0Producer** (EDProducer) - Reconstructs V0 vertices: Ks -> pi
pi, Phi -> KK, D0 -> K pi, D* -> D0 pi, Lambda -> p pi. Configurable
mass windows, track quality, vertex quality, and displacement
significance cuts.

### Muon Identification

**BmmMuonIdProducer** (EDProducer) - Adds custom muon ID variables
using multiple XGBoost models. For each muon, extracts 28+ features
(pT, eta, track quality, chamber matches) and stores per-model MVA
scores. Also performs L1 muon matching and trigger object matching.

**MuonWithSoftMvaProducer** (EDProducer) - Copies slimmedMuons
collection and recomputes soft MVA ID using
PhysicsTools/PatAlgos/SoftMuonMvaEstimator. Optionally runs POG muon
ID MVA (ONNX model, pT > 10 GeV).

### Generator-Level

**GenBmmProducer** (EDProducer) - Extracts B-hadron decay information
at generator level. Finds B-hadrons in pruned gen particles, matches
final state particles, performs reco-to-gen matching. Outputs "genbmm"
and "gensummary" CompositeCandidate collections.

**GenDstarProducer** (EDProducer) - Generator-level D* -> D0(-> K pi)
pi reconstruction.

### Filters

**BmmProdFilter** (EDFilter) - AOD/MiniAOD level dimuon event filter.
Selects loose, tracker, high-purity muon pairs by distance of closest
approach and kinematic cuts.

**MuonFakeFilter** (EDFilter) - Selects events where muons match
hadrons (pions/kaons/protons) at gen level. For fake rate studies.

**BxFilter** (EDFilter) - Object count filter for CompositeCandidate
collections.

### Monitoring and Utilities

**PrimaryVertexInformation** (EDProducer) - Stores primary vertex
properties: position, errors, fit quality, vertex score, sum pT/pT^2.

**TriggerPrescaleProducer** (EDProducer) - Stores trigger prescale
values for selected HLT paths, matched by regex.

**TrackOccupancyLumiProducer** (EDProducer, global with
stream/lumi caching) - Monitors track occupancy per luminosity block
in eta-phi bins. Accumulates across events, merges streams at lumi
boundaries.

**TrackerSeededPFCandAnalyzer** (EDAnalyzer) - Studies tracker-seeded
PF candidates for tracking efficiency.

**SimpleFlatTableProducerPlugins** - Template instantiations for
NanoAOD FlatTable producers (reco::Candidate,
pat::CompositeCandidate, GenEventInfoProduct, HTXS).

## Python Configuration

### nano_cff.py - Entry Point

Defines the CMSSW customization functions called via `--customise`:

- `nanoAOD_customizeDileptonPlusX()` - Main workflow. Loads
  DileptonPlusX_cff and UpdateSlimmedMuons_cff, adds sequences for
  data/MC, enables full genParticle storage.
- `nanoAOD_customizeV0ForMuonFake()` - V0 reconstruction workflow.
- `nanoAOD_customizeBmmMuonId()` - Muon ID workflow.
- `run3_nanoAOD_124()` - Run 3 era adjustments.

### DileptonPlusX_cff.py

Configures DileptonPlusXProducer with all input collections, kinematic
cuts (muon pT > 2, kaon pT > 0.5, pion pT > 0.3), mass windows,
displacement significance thresholds (sigLxy > 3), MVA model paths,
BDT training files. Defines MC variant and NanoAOD table definitions.

Pre-selection knobs for the LLX (llk, llkk) loop:
- `minBhhllVtxProb` (default 0.001): skip B -> ll+hh candidates whose
  ultimate 4-track vtx probability is below threshold.
- `minLLSigLxyForBLLX` (default -1, disabled): require the dilepton
  SV Lxy-significance before spending time on the X-side loop.
- `minHadIPSigBSForBLLX` (default -1, disabled): require the hadron
  track IP/sigma with respect to the beam spot.

Reconstruction modes are individually toggleable: `recoElElX`,
`recoElMu`, `recoMuMuGamma`, `recoDstar`, `recoD0pipi`, `recoD0Kpi`,
`recoKspipi`, `recoKstar`, `recoJpsiKsSlimmed` (B -> J/psi Ks via
slimmedKshortVertices).

Per-table precision is controlled by three named constants used in
every `Var(..., precision=...)` call:
- `full_precision = -1` - stored as full float32 (~7 significant
  digits). Used for masses, pt/eta, vtx_prob, lxy, sigLxy, l3d,
  sl3d, all alphas/cosAlphas.
- `medium_precision = 12` - ~3.5 significant digits. Used for
  vertex errors, isolation variables, impact parameters, lifetime,
  sumpt, most secondary quantities.
- `low_precision = 6` - ~1.5 significant digits. Used for track
  quality variables and MVA scores.

The same constants are defined and applied in `BmmV0ForMuonFake_cff.py`
and `BmmMuonId_cff.py`. File-size reduction comes from the
compression of truncated mantissas; physics-critical observables stay
at full precision.

### BmmMuonId_cff.py

Configures BmmMuonIdProducer with XGBoost model paths (from
data/muon_mva/), 28+ feature definitions, 20+ HLT trigger names, L1
matching. Has era modifiers for Run 2 vs Run 3 trigger collection
names.

### BmmV0ForMuonFake_cff.py

Configures BmmV0Producer with mass windows for Ks, Phi, D0, D*,
Lambda, Ds. Defines output FlatTables for each decay mode.

### Other Configs

- `UpdateSlimmedMuons_cff.py` - MuonWithSoftMvaProducer config
- `ScoutingDileptonPlusX_cff.py` - Scouting workflow config
- `MuonFakeFilter_cfi.py` - Filter parameters (minPt=4, maxEta=1.4)
- `BxToMuMuFilter_cff.py` - Count filter (minNumber >= 1)
- `triggerFilter_cfi.py` - HLT path filter
- `selection.py` - ROOT-to-Python selection expression converter

## Data Files

### ML Models (data/)

**TMVA BDT** - Three event-split models for dilepton classification:
`Run2017-2018-20200515-*-Event{0,1,2}.model` (~27 MB each) with
associated `.features` and `.params` files.

**XGBoost** (data/muon_mva/) - Soft muon MVA models, e.g.
`Run2022-20231030-1731-Event0.model` with `.features` file listing
input variables.

**ONNX** - `mvaID.onnx` (3.7 MB) for POG muon ID MVA.

### Certification (data/certification/)

Golden JSON files for data quality filtering.

## Data Flow

```
MiniAOD / Scouting Input
         |
         v
  UpdateSlimmedMuons          (recompute soft MVA)
         |
         v
  BmmMuonIdProducer           (add XGBoost muon MVA)
         |
         v
  DileptonPlusXProducer       (main dilepton reconstruction)
  BmmV0Producer               (V0 vertex reconstruction)
  PrimaryVertexInformation    (PV properties)
  GenBmmProducer              (MC gen-level info)
  TriggerPrescaleProducer     (prescale storage)
  TrackOccupancyLumiProducer  (occupancy monitoring)
         |
         v
  NanoAOD FlatTables -> ROOT output
```

## Key Algorithms

### Kinematic Vertex Fitting

1. Build transient tracks from reco::Track
2. Create KinematicParticles with mass hypotheses
3. Fit vertex with KinematicParticleVertexFitter or
   KinematicConstrainedVertexFitter (with TwoTrackMassKinematicConstraint
   for J/psi- or phi-constrained fits)
4. Optionally apply mass/pointing constraints
5. Extract mass, position, momentum, chi2/ndof
6. Post-process for Lxy, alpha, decay time with error propagation
7. `set_tree()` guards against NaN/Inf state via `hasFiniteValues()`;
   ill-conditioned fits that "converge" with garbage are invalidated.

### Pre-fit gates (fitBToKLL and fitBToLLhh)

The direct 3-track (mumu+K) and 4-track (mumu+hh) fits are sensitive to
near-degenerate configurations; to reject these upfront and avoid
heap-corrupting matrix inversions, both functions first run a cheaper
composite-particle fit as a stability gate:

fitBToKLL (mumu+K):
1. Run 2-track mumu vertex fit -> composite J/psi KinematicParticle.
2. If requested, apply `MassKinematicConstraint(JPsiMass, JPsiMassErr)`
   to the composite via KinematicParticleFitter.
3. Combine composite J/psi + kaon TransientTrack -> 2-particle fit.
4. Require `valid()` and `vtxProb() > 0.001`, else return empty.
5. On pass, run the direct 3-track ultimate fit (unchanged).

fitBToLLhh (mumu+hh):
1. Run 2-track mumu vertex fit -> composite J/psi.
2. Run 2-track hh vertex fit -> composite hh.
3. Apply `MassKinematicConstraint` to either composite if
   ll_mass_constraint or hh_mass_constraint is requested.
4. Combine the two composites -> 2-particle fit.
5. Require `valid()` and `vtxProb() > 0.001`, else return empty.
6. On pass, run the direct 4-track ultimate fit (unchanged).

The gate cuts combinatorial (non-mc-matched) background by ~50-60%
and preserves tight signal selections; the loose-level mc-matched
signal drops a few percent. Crashes tied to ill-conditioned fits no
longer reproduce.

### B -> J/psi Ks using slimmedKshortVertices

When `recoJpsiKsSlimmed=True`, `buildBToJpsiKsCandidates()` iterates
`slimmedKshortVertices` rather than the PF-track combinatorial loop.
For each Ks vertex it extracts the two daughter PackedCandidates,
builds a local FitCandidate (pion+pion hh vertex), then invokes
`fitBToJpsiKs()`: mumu + composite Ks -> 3-particle fit with optional
J/psi and Ks mass constraints. Output: `BToJpsiKs` collection. This
bypasses the gen-level track-matching problem of Ks daughters whose
tracks start far from the primary vertex.

### Ks selection tiers

`BuildKsCandidate(...)` is the raw factory returning a FitCandidate or
nullptr. It supports three selection tiers via `KsSelectionType`:
- `LooseKsSelection` - mass in [minKsMass, maxKsMass] only.
- `NominalKsSelection` - `ks_loose` cuts plus track IP/sigma > 1.
- `KsmmSelection` - `ks_loose` or `ks_sideband`
  (|lxy|>1 cm, wider `minKsMassVeryLoose/maxKsMassVeryLoose` window),
  plus minimum track pt and eta.

`buildAndDressKsCandidate(...)` is the dressed-output wrapper that
also runs `fillIsolationInfo` and appends to the output collection.

### Dilepton Reconstruction Pipeline

1. Select quality muons/electrons (pT, eta, track quality)
2. Form all dilepton combinations, fit vertices
3. Add optional hadrons/photons for extended modes
4. Compute displacement (Lxy significance, alpha angle)
5. Calculate isolation (close track counts at various thresholds)
6. Evaluate BDT/XGBoost discriminants
7. Gen-match for MC (dR + dPt criteria)
8. Store as pat::CompositeCandidate with user variables
9. Convert to NanoAOD FlatTable format

### Physics Quantities Computed

**Vertex**: position, covariance, chi2, ndof, probability, sum
pT/pT^2

**Displacement**: Lxy (error, significance), alpha angle, distance
of closest approach, longitudinal IP, decay time

**Kinematics**: invariant mass (error), refitted daughter momenta,
p3/p4

**Isolation**: close track counts by vertex probability, displacement
significance, and SV/PV match quality; minimum DOCA

**MVA scores**: TMVA BDT discriminant, XGBoost soft muon ID, POG
muon ID

## Postprocessing System

Located in `postprocess/`. Produces flat ntuples from NanoAOD output
for final physics analysis.

**Processor classes** (inherit from Processor base in
PostProcessingBase.py):
- `Skimmer` - nanoAOD-tools based event filtering
- `SimpleSkimmer` - RDataFrame-based with lumi mask support
- `FlatNtupleFor*` - 8+ specialized flat ntuple producers for
  different analysis tasks (BmmMva, MLFit, MuonMVA, TrigEfficiency,
  DstarFit, Ksmm, etc.)

**Job management**: JobCreator scans EOS input, groups files by
dataset, generates JSON job configs. JobDispatcher schedules across
LocalResourceHandler (vocms118, 16 slots) and SSHResourceHandler
(~13 remote machines, 24-64 slots each, ~350 total slots).

## Validation Suite

Located in `validation/`. 86+ standalone Python/ROOT scripts covering:

- Signal reconstruction efficiency (gen-level and reco-level)
- HLT and L1 trigger efficiency and turn-on curves
- MVA ROC curves and signal/background separation
- Background composition studies per decay mode
- Muon fake rate analysis by particle source
- Mass resolution and yield extraction
- Data/MC comparison
- Pileup dependency studies

Uses both legacy ROOT and modern RDataFrame approaches.
