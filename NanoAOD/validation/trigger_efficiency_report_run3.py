""" Formatted trigger efficiency report

Group efficiency results and compute their ratios to
be used in the analysis note. All results are extractedd
from result/summary.json

"""

import json
from math import *
output = "results/summary.json"
results = json.load(open(output))

hlt = "HLT_DoubleMu4_3_LowMass wrt HLT_Mu4_L1DoubleMu - "
l1_1 = "HLT_Mu0_L1DoubleMu wrt HLT_Mu3_PFJet40 - "
l1_2 = "HLT_Mu0_L1DoubleMu wrt HLT_Mu8 - "

table = [
    ( hlt,         (hlt  + "Run2022", hlt  + "2022 MC")),
    ( l1_1,        (l1_1 + "Run2022", l1_1 + "2022 MC")),
    ( l1_2,        (l1_2 + "Run2022", l1_2 + "2022 MC")),
    ( hlt,         (hlt  + "Run2023", hlt  + "2022 MC")),
    ( l1_1,        (l1_1 + "Run2023", l1_1 + "2022 MC")),
    ( l1_2,        (l1_2 + "Run2023", l1_2 + "2022 MC")),
] 

def ratio(data, mc):
    r = data['eff'] / mc['eff']
    err = sqrt(pow(data['eff_err'] / data['eff'], 2) + pow(mc['eff_err'] / mc['eff'], 2))
    return (r, r * err)

def print_results(data, mc):
    if data in results:
        print("& $%5.2f \pm %4.2f$ " % (results[data]['eff']*100, results[data]['eff_err']*100), end=' ')
    else:
        print("&                  ", end=' ')

    if mc in results:
        print("& $%5.2f \pm %4.2f$ " % (results[mc  ]['eff']*100, results[mc  ]['eff_err']*100), end=' ')
    else:
        print("&                  ", end=' ')

    if data in results and mc in results:
        print("& $%5.3f \pm %5.3f$ " % ratio(results[data], results[mc  ]), end=' ')
    else:
        print("&                  ", end=' ')


for ( row, (data, mc) ) in table:
    print("%- 30s " % (row), end=' ')
    print_results(data, mc)
    print(r'\\')
