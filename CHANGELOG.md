# Changelog

## New developments

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
