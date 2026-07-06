# Pixel auto-masking in B-physics ntuples

Design and reference for the pixel-quality information added to the Bmm NanoAOD: per-event pixel
ROC quality counts (`pix_*`) from the `SiPixelQuality` conditions, and a per-muon pixel-crossing
status (`MuonId_pixel*`) from the track `HitPattern`. This document explains what the branches
mean, how they are computed, and how to interpret them. Provenance markers: `[verified]` = checked
against `cms-sw/cmssw` sources (list in §6; APIs stable for ≥ CMSSW_14_0_X).

---

## 1. Motivation and background

### 1.1 What auto-masking is

TBM cores intermittently enter a stuck state during Run 3 (typically attributed to SEUs; rate
grows with instantaneous luminosity and radiation damage). The FED detects the dead channel and
auto-masks it until a reset recovers it; in raw data this is FED error type 25 —
`"Error: Disabled FED channel (ROC=25)"` `[verified: SiPixelRawDataError.cc]`. Masked channels
produce no hits. BPix L1 (r ≈ 2.9 cm, PROC600 chips, highest rate/fluence) is the most affected
region; losing L1 hits degrades IP resolution, vertexing, and soft-track efficiency. Masking is
recorded at **ROC granularity** (whole readout chips or groups thereof — never individual pixels).

Two views exist in CMSSW and this design uses both: **per event** (the unpacker records each
event's disabled FED channels, folded into tracking — §5.5) and **per lumisection** (the PCL
conditions payloads — below).

### 1.2 The PCL SiPixelQuality workflow `[verified: SiPixelStatusProducer/Manager/Harvester]`

Input: `ALCARECOSiPixelCalZeroBias`. Per LS, `SiPixelStatusProducer` accumulates per-ROC digi
occupancy and builds the FEDerror25 channel list as the **intersection over all events in the LS**
(a channel is flagged only if off in *every* ZeroBias event of that LS); one status object per LS.
`SiPixelStatusHarvester` (end of run) builds and uploads:

| flavor | content (from the source) | central tag |
|---|---|---|
| `stuckTBM` | FEDerror25 ROCs **minus** permanently-bad ROCs | `SiPixelQuality_byPCL_stuckTBM_v1` |
| `prompt` | permanent bad + "other"; FED-25 ROCs excluded by construction | `SiPixelQuality_byPCL_prompt_v2` |
| `other` | low-occupancy ROCs minus FED25 minus permanent | `SiPixelQuality_byPCL_other_v1` |

The `stuckTBM` flavor keeps **1-LS granularity** (IOV written at the exact LS where the set
changes; lists replaced, not accumulated — recoveries drop out). A ROC that is both permanently
dead and stuck appears only in `prompt`, so the stuckTBM count is the *incremental* loss on top of
the permanent state. `stuckTBM` and `prompt` are **disjoint at ROC level** on Prompt tags.

**Tag facts (conddb, 2026-07-02):** `stuckTBM_v1` is the only production stuckTBM tag and is in
**no GlobalTag** — it must be loaded explicitly (§3.3). It covers runs **315690 → 404926** (all of
2018-late onward; ≈203k IOVs). `stuckTBM_v1` and `prompt_v2` are the live PCL streams, populated in
lock-step; `prompt_v4` is a frozen single-IOV snapshot — do not use. Caveat: `conddb list` silently
caps at the last 500 IOVs unless `--limit` is raised.

### 1.3 IOV mechanics

`since = (run << 32) | LS`; python `run, ls = divmod(since, 2**32)`. **End-of-run closure
asymmetry** `[verified]`: at `lastLS+1` the harvester closes `prompt` (reverts to permanent-bad)
and `other` (empty) — **not `stuckTBM`**, whose last payload persists into later
non-PCL-processed (run,LS) ranges. Interpret stuckTBM values only for PCL-processed pp runs.

### 1.4 How representative is LS-level information? (1 LS ≈ 23 s)

Because of the full-LS intersection criterion the failure modes are enumerable and **one-sided**
(all under-counting): (i) the *onset* LS is not flagged; (ii) the *recovery* LS is not flagged;
(iii) episodes shorter than ~1–2 LS are invisible. Events labelled "masked" are pure; contamination
sits in the "good" sample. Magnitude, from the tag alone:
`fraction of mislabelled LS ≈ (# stuckTBM IOV boundaries in run) / (# LS in run)`. Two systematics:
onsets cluster at high instantaneous luminosity (don't let the masked count double as a PU proxy),
and the invisible short-episode class needs a per-event cross-check to bound. **Per-muon physics is
unaffected**: the per-muon flag (§2.3) is event-exact.

**Measured magnitudes** (full-run LS scans of the tag, cross-checked against
`plot_SiPixelBPixQualityMap`):

| run | era / activity | L1 mean masked (fraction) | L2 | L3 | L4 | whole-det mean |
|---|---|---|---|---|---|---|
| 385620 | 2024, typical | 62 ROCs (4.0 %) | 105 (2.9 %) | 65 (1.2 %) | 46 (0.6 %) | 378 |
| 398803 | 2025, median | 262 ROCs (**17.1 %**) | 158 (4.4 %) | 165 (2.9 %) | 33 (0.4 %) | 752 |
| 398027 | 2025, busiest | 399 ROCs (**26 %**, max 35 %) | 166 (4.6 %) | 151 (2.7 %) | 54 (0.7 %) | 1126 |

L1 dominates as a *fraction* in every era, overwhelmingly in 2025. Note: **run boundaries reset the
state** — stuck channels are cleared at run start and masking accumulates through the run, so early
lumisections (LS ≲ 10) are near-clean and unrepresentative; and masked > 0 in essentially 100 % of
in-run LS for L1–L4, so the §2.4 estimators matter everywhere, not just in special runs.

---

## 2. The ntuple branches

### 2.1 Region indexing

`0..3` = BPix L1..L4; `4..6` = FPix minus D1..D3; `7..9` = FPix plus D1..D3 (a fixed 10-region
superset covering both pixel geometries; §3.2).

### 2.2 Event level — the `pix` table

A 10-row `nanoaod::FlatTable` named `pix` (counter `npix` ≡ 10), one row per region, constant
within a conditions IOV. `run`/`luminosityBlock` come from standard NanoAOD.

| branch | type | definition |
|---|---|---|
| `pix_region` | int | region index per §2.1 |
| `pix_nRocsTotal` | int | total ROCs in region from geometry; **0 = region absent in this geometry** (Phase-0 2016 lacks 3/6/9); −1 = not read |
| `pix_nRocsMasked` | int | ROCs in the **stuckTBM** payload at this (run,LS) — auto-masked, permanent part subtracted; −1 if stuckTBM not read |
| `pix_nRocsDead` | int | ROCs in the **base** quality (label `""`; Prompt GTs: permanent + persistently low-occupancy) |
| `pix_nRocsBadTotal` | int | ROC-level **union** of masked and dead (robust for ReReco tags that merge sources); equals dead when stuckTBM not read |
| `pix_nRocsActive` | int | `nRocsTotal − nRocsBadTotal` |
| `pix_nModulesMasked` | int | modules with **all** their ROCs auto-masked; −1 if stuckTBM not read |
| `pix_fRocsMasked` | float | `nRocsMasked / nRocsTotal` (geometric proxy for binning, not an efficiency); −1 if undefined |

**Sentinel convention (three distinguishable states)**: `nRocsTotal == −1` → not read (MC, or
`readPixelQuality=False`); `nRocsTotal == 0` → region absent in the active geometry (Phase-0
regions 3/6/9); `nRocsTotal > 0` → real. The masked columns carry `−1` when the stuckTBM tag is not
read (all Run-2, §3.3) — distinct from a real `0`.

### 2.3 Per-muon level — `MuonId_pixel*`

Added in `plugins/BmmMuonIdProducer.cc` (next to the `fill_track_info` call), declared in
`BmmMuonIdVariables` (`python/BmmMuonId_cff.py`); rows are index-aligned with the standard `Muon`
table, so candidate-level access is `Take(MuonId_x, mm_muN_index)` as other consumers already do
(e.g. `validation/rdf_kmm_trigger_efficiency.py`). Derived from the embedded inner-track
`HitPattern`, so these need **no conditions** and work on every era and data tier.

| branch | type | definition |
|---|---|---|
| `MuonId_pixelL1Status` | int | 0 valid L1 hit; 1 crossed **INACTIVE** L1 area in *this event*; 2 crossed active area, no hit; 3 bad hit; 4 not crossed / no track |
| `MuonId_pixelStatusWord` | int | the same status packed 3 bits/region × 10 regions. **No FPix z-side in HitPattern (§2.5 C-g)** → all FPix crossings appear in regions 4–6; regions 7–9 always read 4 |
| `MuonId_pixelBarrelOffLayers` | int | `hitPattern().pixelBarrelLayersTotallyOffOrBad(TRACK_HITS)` |

The status loops **both** `TRACK_HITS` and `MISSING_INNER_HITS` (an L1 crossing of a track whose
innermost valid hit is on L2+ lands in `MISSING_INNER_HITS`), and distinguishes INACTIVE from
MISSING (precedence valid > inactive > bad > missing). Track fallback: `innerTrack()` →
`bestTrack()` → status 4. The main analysis handle is `MuonId_pixelL1Status == 1` — muon crossed an
auto-masked/inactive L1 region in that exact event — used with `pix_nRocsMasked[0]` for the
per-LS view. (A per-daughter copy inside the `mm` table can be added later by extending
`bmm::fill_track_info`; not done in this iteration since `MuonId` is index-aligned.)

### 2.4 Hit-efficiency estimators (offline, from the branches)

```
eff_std   = N(valid) / (N(valid)+N(missing))                 # reco convention: masked regions excluded
eff_incl  = N(valid) / (N(valid)+N(missing)+N(inactive))     # includes masking losses
loss_mask = N(inactive) / (N(valid)+N(missing)+N(inactive))  # masking-induced hit-inefficiency
```

Attribution of `loss_mask` to auto-masking: compare `pix_nRocsMasked > 0` vs `== 0` lumisections at
fixed `pix_nRocsDead`, or the deferred tier-2 footprint attribution (§7).

### 2.5 Interpretation caveats

* **C-a granularity**: event-level branches are per-(run,LS) step functions (one-sided
  under-count; §1.4). The per-muon flag does **not** have this limitation (event-exact).
* **C-b leakage**: no end-of-run closure for stuckTBM → trust masked counts only for
  PCL-processed pp runs.
* **C-c decomposition**: masked ∩ dead = ∅ on Prompt tags; the `nRocsBadTotal` union protects
  against ReReco tags that merge sources.
* **C-d INACTIVE is cause-agnostic** (automask, dead, HV-off, whole-event pixel-off protection).
* **C-e survivor bias**: per-muon flags exist only for *reconstructed* muons; masking also kills
  seeds/tracks outright — quantify with tag-and-probe or MC quality scenarios. State in the note.
* **C-f hadrons on MiniAOD**: no per-hadron HitPattern; only the deferred tier-2 (§7) could flag them.
* **C-g FPix side** `[verified: HitPattern.cc:178 side=isStereo, :877–897 isStereo==0 for pixels]`:
  the HitPattern "side" bit is the strip **stereo** flag and is 0 for every pixel hit — HitPattern
  carries **no FPix z-side in Run 3**. So in `pixelStatusWord` all FPix crossings land in the
  side-inclusive regions 4–6 and regions 7–9 always read 4. For per-side FPix, use the track η sign
  or tier-2.

---

## 3. How it is computed

### 3.1 Files and wiring

* `plugins/PixelQualityTableProducer.cc` — new `edm::stream::EDProducer<>` emitting the `pix`
  table. ESWatcher-gated: it recomputes counts only when the `SiPixelQualityFromDbRcd` IOV changes,
  otherwise just refills from the cached arrays (modeled on `TriggerPrescaleProducer.cc`).
* `python/PixelQuality_cff.py` — `pixelQualityTable` (+ `pixelQualityMcTable = clone(readPixelQuality=False)`).
* `python/nano_cff.py` — loads the cff, appends the tables to the data / MC sequences, decides the
  stuckTBM tag automatically (§3.3), and carries the Run-2 tau bypass (§3.4).
* `src/CommonTools.cc` / `interface/CommonTools.h` — `bmm::pixel_region_cross_status` and
  `get_pixel_status_word` (the per-muon HitPattern logic).
* `plugins/BmmMuonIdProducer.cc`, `python/BmmMuonId_cff.py` — the `MuonId_pixel*` userInts + Vars.

FlatTable module labels must end in `Table` to match the default NanoAOD keep patterns; the MC twin
tables/sequences are wired in parallel. Per-event cost is negligible (measured: `pixelQualityTable`
≈ 0.015 ms/event, the per-muon logic adds no measurable time — ~0.004 % of a 354 ms/event Bmm job,
well below run-to-run noise; peak memory unchanged; output +26 B/event).

### 3.2 Geometry: Phase-0 vs Phase-1

ROC counts per module come from `TrackerGeometry` / `PixelTopology` at runtime, **not** hardcoded,
so one producer handles both pixel detectors:

* **Phase-1 (2017+, Run 3)**: all modules 16 ROCs. BPix L1–L4 = 96/224/352/512 modules →
  **1536/3584/5632/8192** ROCs; FPix 112 modules/disk/side → **1792** ROCs, ×6 disk-sides.
* **Phase-0 (2016)**: 3 BPix layers + 2 FPix disks/side; modules are **not** all 16 ROCs — BPix
  half-modules and FPix plaquettes give per-module ROC counts {2,5,6,8,10,16}, so BPix L1–L3 totals
  are **2304/3840/5376** and each FPix disk-side is **1080**. Regions 3/6/9 are absent (`nRocsTotal
  == 0`).

**Every ROC has the same 80×52 = 4160 channels** in both detectors (incl. the Phase-1 L1 PROC600),
so ROC counting is a uniform unit and `fRocsMasked` is a fair area proxy. The bad-ROC loop is
bounded to each module's real ROC count, so a wholly-bad sub-16 module (`errorType==0`) is not
over-counted.

### 3.3 Conditions loading (automatic, GT-driven)

The `stuckTBM` tag is in no GlobalTag, so it is loaded via a `GlobalTag.toGet` append under an
additive ES label `"stuckTBM"` on record `SiPixelQualityFromDbRcd`; the base/dead payload is
whatever the job GT already serves at label `""` (Prompt GTs: `byPCL_prompt_v2`). The producer reads
both with one `edm::ESWatcher<SiPixelQualityFromDbRcd>`.

Whether to load stuckTBM is decided **automatically from the job's GlobalTag** —
`pixel_stuckTBM_enabled(gt)` returns `'dataRun3' in gt`; the customise sets
`pixelQualityTable.readStuckTBM` and appends the `toGet` accordingly. Rationale: the tag's first IOV
is run 315690, and Run-3 (runs ≥ 355k) is the only class guaranteed to be covered; **all Run-2**
(2016/2017/2018) and MC → tag skipped, masked columns = −1 sentinel, while dead/total/geometry and
the per-muon flags still work.

The append **must be data-only and coverage-guaranteed**
`[verified: CondDBESSource.cc:567/576–577/598]`: `CondDBESSource::setIntervalFor` **intersects IOV
validity across all labels of a record**, so an appended-but-uncovered `stuckTBM` label would
invalidate the whole `SiPixelQualityFromDbRcd` record (killing the base/dead read too) — fatal.
Hence the conservative, config-time GT gate (never appended for MC — which runs at run ≈ 1 — or for
Run-2). Verified with `edmConfigDump`: the tag appears exactly once on Run-3 data, zero on MC/Run-2.

### 3.4 Run-2 tau-reprocessing bypass

Independent of pixels, running the Bmm NanoAOD on Run-2 (era `run2_nanoAOD_106Xv2`) aborts in the
**standard** NanoAOD tau sequence: it re-embeds the legacy tau anti-electron MVA6 (and MVArun2
isolation MVAs) whose `GBRForest` payloads no longer exist in the conditions DB for ≥ 14_0_X. Bmm
uses no taus, but `PATObjectCrossLinker` (`linkedObjects`, consumed by `BmmMuonId`) needs valid
`pat::Tau` collections. `disable_legacy_tau_reprocessing(process)` (in `nano_cff.py`) points the
final tau selectors at the raw MiniAOD taus (walking each reprocessing producer's `src` back to
source, which also handles the 15_0_X `@skipCurrentProcess` same-name embedder) with a trivial `pt`
cut, and drops the tau tables. It is gated on the presence of an anti-electron-rejection producer,
so it is a **no-op for Run 3**, and it keys on module contents rather than names, so it works in
14_0_X and 15_0_X alike.

It **must be invoked from `--customise_commands`** (which run last), not as a `--customise`
function: cmsDriver builds the tau reprocessing in `nanoAOD_customizeCommon`, which runs *after* the
user `--customise` functions, so a `--customise` call would run too early and be clobbered. The
ready-to-run Run-2 cmsDriver commands (data and MC) are in `Bmm5/README.md` (Processing examples).

---

## 4. Verified CMSSW reference

### 4.1 Payload class `SiPixelQuality` (`CondFormats/SiPixelObjects`)

```cpp
struct disabledModuleType {
  uint32_t DetID;          // offline DetId of the module
  int errorType;           // 0="whole" module bad, 1="tbmA", 2="tbmB", 3="none" (use BadRocs)
  unsigned short BadRocs;  // 16-bit mask, bit n set <=> ROC n bad
};
const std::vector<disabledModuleType> getBadComponentList() const;  // returns BY VALUE (copy)
```

errorType ↔ BadRocs: 0/whole→65535, 1/tbmA→255, 2/tbmB→65280, 3/none→arbitrary mask
`[verified: SiPixelQualityESProducer.cc]`; the PCL sets errorType=0 with the full mask.
**Trap** `[verified: SiPixelQuality.cc]`: `IsModuleBad` / `IsRocBad` / `getBadRocs` each **copy and
re-sort the whole payload per call** — never use them in per-ROC loops; iterate
`getBadComponentList()` once per IOV.

### 4.2 Records

* `SiPixelQualityFromDbRcd` — raw DB payload record; what the PCL tags populate; **read this**.
* `SiPixelQualityRcd` — dependent record merging the DB payload with `SiPixelDetVOffRcd` HV-off
  info; used by reco; do **not** use for the masked/dead decomposition.

### 4.3 Topology / geometry helpers

`TrackerTopology` (`TrackerTopologyRcd`): `pxbLayer` 1–4; `pxfSide` 1=minus/2=plus, `pxfDisk` 1–3.
`PixelSubdetector::PixelBarrel=1, PixelEndcap=2`. `TrackerGeometry` (`TrackerDigiGeometryRecord`):
`detsPXB()`, `detsPXF()`; `PixelGeomDetUnit::specificTopology().rocsX()*rocsY()` = ROCs/module,
`rowsperroc()*colsperroc()` = channels/ROC.

### 4.4 The per-event chain: auto-masking → track HitPattern `[verified]`

1. Unpacker produces per event `PixelFEDChannelCollection` (label `siPixelDigis`; element
   `{fed, link, roc_first, roc_last}`) listing channels off in *that event*, incl. auto-masked.
2. `MeasurementTrackerEventProducer` wires it into tracking
   (`badPixelFEDChannelCollectionLabels = ['siPixelDigis']`).
3. `TkPixelMeasurementDet::hasBadComponents()` tests the predicted local position (3σ window)
   against SiPixelQuality bad-ROC rectangles **and** the per-event bad-FED spans; if inside, the
   crossing is recorded **INACTIVE** instead of MISSING.
4. `reco::HitPattern` stores per crossing: subdet, layer/disk, side, hit type
   `VALID=0, MISSING=1, INACTIVE=2, BAD=3`, in categories `TRACK_HITS` / `MISSING_INNER_HITS` /
   `MISSING_OUTER_HITS`. Key API: `numberOfAllHits(cat)`, `getHitPattern(cat,i)`,
   `pixelBarrelHitFilter`/`pixelEndcapHitFilter`, `getLayer`, `getHitType`,
   `pixelBarrelLayersTotallyOffOrBad(cat)`. `numberOfLostHits` counts only MISSING (excludes
   INACTIVE) — hence the explicit hit-type loop.
5. MiniAOD keeps it: `patMuons.embedTrack = True`, and `HitPattern` is a data member of
   `reco::Track` → `mu.innerTrack()->hitPattern()` works on MiniAOD.

### 4.5 Data-tier availability

| object | RAW | RECO | AOD | MiniAOD |
|---|---|---|---|---|
| `PixelFEDChannelCollection` (`siPixelDigis`) | re-unpack | kept | no | no |
| muon inner-track `HitPattern` (incl. INACTIVE) | — | yes | yes | **yes** (embedded) |
| hadron-track full `HitPattern` | — | yes | yes | **no** (`lostInnerHits()` summary only) |
| `SiPixelQuality` conditions | any tier (EventSetup) | | | |

### 4.6 conddb cheat sheet

```bash
conddb list SiPixelQuality_byPCL_stuckTBM_v1 --limit 1000000 | tail   # IOV coverage (run:LS)
conddb listGTsForTag SiPixelQuality_byPCL_stuckTBM_v1                  # which GTs contain it
conddb list <GT> | grep -i SiPixelQuality                             # what the job GT gives at label ""
```

For per-region validation against the payload inspector, use
`plot_SiPixelQualityBadRocsTimeHistory` (whole detector) or the `plot_SiPixel{BPix,FPix}QualityMap`
maps. NB: the per-layer `_L1.._D3` TimeHistory variants do **not** exist in CMSSW_14_0_16 (master
only); use the whole-detector plot or the producer's `verbosePixelQuality` log.

---

## 5. Troubleshooting

* `No data of type SiPixelQuality ... label "stuckTBM"` → the `toGet` append is missing (the GT was
  not recognized as Run-3 data; see §3.3).
* `Fatal Exception ... Forward compatibility cannot be supported` opening 2025 data → the file was
  written with CMSSW_15_0_X; use a ≥ 15_0_15 area (the package builds there unchanged).
* `No data of type GBRForest ... RecoTauTag_antiElectronMVA_...` on Run-2 → the tau bypass
  (§3.4) was not applied; it must be in `--customise_commands`, not `--customise`.
* Masked counts ≈ 0 in the first few LS of a run → real: the stuck-TBM state resets at run start and
  accumulates through the run (§1.4), not a bug.
* Nonzero masked counts in a non-collision run → end-of-run leakage (§1.3 / C-b), not a bug.
* Muon status always 4 → track refs not embedded, or the module runs before pat muons exist.
* Checking whether the `toGet` reached a job: grep the **`edmConfigDump` output** for the tag name
  `SiPixelQuality_byPCL_stuckTBM_v1` (expect 1 on Run-3 data, 0 on MC/Run-2).

---

## 6. Appendix — verified sources (cms-sw/cmssw)

* `CondFormats/SiPixelObjects/interface/SiPixelQuality.{h,cc}` — payload structure; errorType
  semantics; accessor copy+sort trap.
* `CalibTracker/SiPixelESProducers/plugins/SiPixelQualityESProducer.cc` — errorType↔BadRocs table.
* `CalibTracker/SiPixelQuality/plugins/SiPixelStatusProducer.cc`, `SiPixelStatusHarvester.cc`,
  `src/SiPixelModuleStatus.cc` — full-LS intersection; stuckTBM = FEDerror25 − permanent;
  IOV-on-change; end-of-run closure for prompt/other only.
* `CondCore/SiPixelPlugins/plugins/SiPixelQuality_PayloadInspector.cc` — reference counting; plots.
* `CondCore/ESSources/plugins/CondDBESSource.cc` — `toGet` label support; cross-label IOV
  intersection (§3.3).
* `DataFormats/SiPixelDetId/interface/PixelFEDChannel.h`; `DataFormats/SiPixelRawData/src/SiPixelRawDataError.cc`.
* `RecoTracker/MeasurementDet/plugins/{MeasurementTrackerEventProducer,TkPixelMeasurementDet}.cc` —
  per-event FED channels wired into tracking; INACTIVE classification.
* `DataFormats/TrackReco/interface/HitPattern.h`, `src/HitPattern.cc` — hit types/filters; side bit
  = strip stereo (§2.5 C-g).
* `PhysicsTools/PatAlgos/python/slimming/miniAOD_tools.py` — `patMuons.embedTrack = True`.
* `Geometry/CommonTopologies/interface/PixelGeomDetUnit.h`, `PixelTopology.h` — ROCs/module,
  channels/ROC.

---

## 7. Remaining work

* **Tier-2 masked-footprint attribution** (deferred): per-track propagation to masked-ROC
  rectangles (`SiPixelFedCablingMap` + propagator), giving a cause-resolved `MuonId_pixelL1MaskedCross`
  and the only per-candidate option for B-decay **hadrons** on MiniAOD.
* **Full-statistics V3**: the masked-fraction ↔ `MuonId_pixelL1Status==1` correlation on
  production-scale ntuples (initial pass done: L1-inactive rate 3.6 % / 15.3 % / 25.2 % at
  0 % / ~11 % / 31 % masked L1).
* **2025 data GT**: add a 2025 prompt GT + a ≥ 15_0_15 test sample to `test/test_config.py`
  (the repo currently pins only the 2024 data GT).
