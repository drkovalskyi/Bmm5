# Changelog

## [NanoAODv15-V01] - 2025/10/28
### Changes
	- added data certification (muon quality)
	- disabled warning from TwoTrackMinimumDistance
	- added support for BtoLLhh instead of just BtoLLkk
	- added kinematic fits for TnP

## [NanoAODv14-V06] - 2025/09/20
### Changes
	- Implemented track reco efficiency measurement with stand alone muons
	- Added STA muons to the selection and store extra information for them
	- Added inner track algorithm mask

## [NanoAODv14-V05] - 2025/09/02
### Changes
	- Added loose Ks to two charged pf candidates with lxy>1.

## [NanoAODv14-V04] - 2025/04/24
### Changes
	- Lowered minimum muon pt to 2 GeV for the muon collection

## [NanoAODv14-V03] - 2025/01/10
### Changes
	- added composite candidate gen matching for K_S

## [NanoAODv14-V02] - 2024/12/29
### Changes
	- save tracking information for dileptons and hh final states
	- relaxed Ks mass preselection and applied main selection on refitted mass
	- added BsToPhiPhi

## [NanoAODv14-V01] - 2024/10/11
### Changes
	- add displaced tracks to classify dimuon PV
	- enabled tracking isolation for Kspipi to use as a control sample for Ksmm
	- fixed L1 matching
	- switched to CMSSW_14_0_X
	- added muon PV association
	
## [NanoAODv12-V06] - 2024/05/21
### Changes
	- bug fix for invalid pv index

## [NanoAODv12-V05] - 2024/05/10
### Changes
	- Added kstar->ks pi
	- Added HLT scouting prototype
	- Added tracks and isolation information
	- Added pixel pattern information for muons and tracks

## [NanoAODv12-V04] - 2023/12/15
### Changes
	- Added tripple muon final state
	- opened up muon fakes phase space for muon id paper (need to reduce it later)

## [NanoAODv12-V03] - 2023/12/12
### Changes
	- Latest run3 muon id
	- added gen index to GenPart for some of the matched
	objects to be able to investigate background decays
	- keep all pruned gen particles from MiniAOD by default
	- lowered minimum hadron pt to 3GeV to match the muon trigger 

## [NanoAODv12-V02] - 2023/11/13
### Changes
	- bug fixes for XGBoost
	- Added BsToDsMuNu

## [NanoAODv12-V01] - 2023/10/18
### Changes
	- Migrated to NanoAODv12 setup to include Run
	- Switched to the new soft mva id
	- Added KsToPiPi similar to KsToMuMu
	- Widen D0mm mass window to include more background

## [NanoAODv10-V05] - 2023/08/18
### Changes
	- lowered minimal muon pt to 2 GeV for all types of muons

## [NanoAODv10-V04] - 2023/07/21
### Changes
	- added extra muon variables
	- added Muon POG mva id trained for muons with Pt > 10 GeV
	- added e-mu final state
	- lowered minimal muon pt to 2 GeV
	- bug fixes

## [NanoAODv10-V03] - 2023/05/04
### Changes
	- added primary vertex information
	- refactored code handling displacement calculations
	- added primary vertex fit excluding signal tracks
	- added PV refit information for Dstar decays
	- bug fixes

## [NanoAODv10-V02] - 2023/03/26
### Changes
	- Revised code organization and handling of different final states
	- Re-wrote Dstar code
	- Renamed a number of branches
	- added more info about L1 muons

## [NanoAODv10-V01] - 2022/12/22
### Fixed
	- Improved multithreading support (stream level)
	- Disabled slow options not relevant for Run3 at the moment

## [NanoAODv9-V04] - 2022/12/09
### New
	- Added Jpsi to mumu reconstructed as tracks, i.e. for muons below 4 GeV

## [NanoAODv9-V03] - 2022/10/14
### New
	- Major code re-organization with re-writes
	- Added dielectron channels (Kee, KKee)
	- Added Dstar to D0(mumu) pi

## [NanoAODv9-V02] - 2022/10/14
### New
	- Reorganized code
	- Integrated photon KinematicParticle code
	- Added BtoMuMuGamma final state using regular and converted photons

### Fixed
	- fixed DstarToD0Pi code
	- updated analysis scripts

## [NanoAODv8-V05] - 2021/12/06
### New
	- Added L1 and HLT objects matched to offline muons
	- Added HLT filter information for matching objects

## [NanoAODv8-V04] - 2021/11/05
### New
	- Improved MVA matching between JpsiK and Bmm

### Fixed
	- Gen particle identification in genbmm block is more robust now


## [NanoAODv8-V03] - 2021/08/26
### Fixed
	- Added protection from getting Nan for kin_alpha

## [NanoAODv8-V02] - 2021/07/03
### Added
	- Trigger prescale for selected set of HLT triggers
	- GenFilterInfo in LuminosityBlocks to compute proper event yield normalization

### Fixed
	- GenInfo extraction for JpsiK
	- Added exception handling for bad fits

## [NanoAODv8-V01] - 2021/05/26
	- Ultra Legacy equivalent of NanoAODv6-V17. Requires CMSSW_10_6_19_patch2

## [NanoAODv6-V17] - 2021/05/10

### Fixed
	- Added protection from getting NaN for pointing angle in a special case

## [NanoAODv6-V16] - 2021/05/07

### Fixed
	- Removed confusing pointing angle definitions and renamed some of them
	- Added protection from NaN in impact parameter calculation

### Added
	- new soft muon MVA
	- uncertainty calculation for pointing angle
	- reorganized/modified some core code - bugs are possible

## [NanoAODv6-V15] - 2021/04/10

### Fixed
	- Fixed bug in JpsiKK code
	- Increase DsToPhiPi mass windows

## [NanoAODv6-V14] - 2021/04/02

### Fixed
	- use only a single thread in XGBooster

### Added
	- muon id inputs for new MVA id
	- DsToPhiPiToKKPi control samples
	- expanded JpsiK to mmK including Psi(2S)K
	- better postprocessing

## [NanoAODv6-V13] - 2021/01/05

### Fixed
	- catch std::exceptions during fitting and mark fit as invalid

## [NanoAODv6-V12] - 2020/12/08

### Fixed
	- bug fix in Lambda reconstruction

## [NanoAODv6-V11] - 2020/12/04

### Added
	- V0 impact parameter info
	- extra impact parameter info for B candidates

### Fixed
	- enforce PF candidate mass hypothesis in V0 candidates

## [NanoAODv6-V10] - 2020/11/27

### Added
	- Added D0ToKPi, LambdaToPPi and PhiToKK for muon fake studies
	- Added track impact parameter for V0

### Modified
	- MC matching for V0 is modified to inmprobe matching efficiency

## [NanoAODv6-V09] - 2020/08/21

### Fixed
	- Fixed bug in MuonCand build from a hadron (the track reference could be wrong)

## [NanoAODv6-V08] - 2020/08/06

### Added
	- Injected hadrons from relevant exclusive decays into good muon candidates to study muon fakes. This may break some old code since muon index for hadron is a negative number and cannot be used to reference a muon unless it was reconstructed as a muon
	- Added gen summary for all relevant exclusive backgrounds
	- Added a few variables
### Modified
	- Redesinged genbmm block. It is not compatible with the old code.
	- Improved gen matching (adde extra information and started using packed candidates)
### Fixed
	- Fixed bug in matching tracks to ignore

## [NanoAODv6-V07] - 2020/05/26

### Added
	- Recompute soft muon MVA

## [NanoAODv6-V06] - 2020/05/15

### Added
	- Bmm emulation with BtoJpsiK events
	- skimming tools
	- KsToPiPi reconstruction for muon fake studies
	- dimuon vertexing with pointing constrain (kinpc)
	- XGBoost is integrated in NanoAOD code to produce ntuple with new MVA centrally
	- Generator level filter for Bmm signature to effectively select relevant QCD events
	- MC production config files

## [NanoAODv6-V05] - 2020/03/16

### Fixed
	- Fixed bug in generator level information extraction critical for efficiency studies

## [NanoAODv6-V04] - 2020/03/06

### Added
	- kk mass for BsToJpsiPhi
	- Bmm gen event level information for efficiency studies
	- event classification based on the production type
	- decay time in 3D and 2D
	- gen decay time
	- added some missing kinematic variables

### Fixed
	- fixed refitted daughter information
