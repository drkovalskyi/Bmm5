# Soft Muon MVA selectors

## Run 3

New Soft Muon MVA model is trained using 800M InclusiveDileptonMinBias
MC simulated events. The muon candidates originating from the pion and
kaon decays in flight are classified as background and thoses matched
to the heavy flavour decays as signal. The matching is performed using
the sim-hit information in the muon detectors. The training sample
includes muon candidates isTrackerMuon or isGlobalMuon requirements
and have pt greated than 2GeV.

The model is trained using muon pt and eta among input features
without event reweighting. This approach can be sub-optimal for some
analyses.


### Models

* Run2022-20231030-1731-Event0
  * eta and pt re-weighted
  * selection: pt>2, TrackerMuon or GlobalMuon
  * additiona inputs: "glbNormChi2", "trkLayers", "highPurity"
* Run2022-20230926-1028-Event0
  * eta and pt re-weighted
  * selection: pt>2, TrackerMuon or GlobalMuon
* Run2022-20231115-1932-Event0
  * eta and pt re-weighted
  * selection
     * pt>2 GeV
     * TrackerMuon or GlobalMuon
     * HLT trigger
* Run2022-20231113-0145-Event0.
  * selection
     * pt>2 GeV
     * TrackerMuon or GlobalMuon

## Run2018

Files starting with Run2018 contain the soft mva id trained with Run2
using the XGBoost BDT algorithm. It performed better than the official
MVA based on Run2016 TMVA BDT, but the gain was not significant enough
to warrant an new review process.