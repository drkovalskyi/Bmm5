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

## Run2018

Files starting with Run2018 contain the soft mva id trained with Run2
using the XGBoost BDT algorithm. It performed better than the official
MVA based on Run2016 TMVA BDT, but the gain was not significant enough
to warrant an new review process.