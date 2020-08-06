# Changelog

## New developments

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
