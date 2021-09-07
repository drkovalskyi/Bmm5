# Training MVA with XGBoost

## Code organization

The code works for a wide variaty of tasks. The core is in
ModelHandler.py, which is a wrapper around XGBoost. All that you need
to do is to create a new class that inherits from ModelHandler,
provide feature names, input files and tune the parameters.

`plot_results_*` are plain Python scripts to make comparison plots

## Prepare Input Datasets

We use flat ntuples produced with
https://github.com/drkovalskyi/Bmm5/blob/master/NanoAOD/postprocess/FlatNtupleForBmmMva.py
code

## Examples

`python train_bmm_mva.py`

`python plot_results_bmm_mva.py`
