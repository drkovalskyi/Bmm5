""" Formatted trigger efficiency report

Group efficiency results and compute their ratios to
be used in the analysis note. All results are extractedd
from result/summary.json

"""

import json
from math import *
output = "results/summary.json"
results = json.load(open(output))

jpsik1 = "HLT_DoubleMu4_3_Jpsi_Displaced wrt HLT_Dimuon0_LowMass_L1_0er1p5 trig eff for Jpsi - "
jpsik2 = "HLT_DoubleMu4_3_Jpsi wrt HLT_Dimuon0_LowMass_L1_0er1p5 displaced trig eff for Jpsi - "
jpsik3 = "HLT_DoubleMu4_3_Jpsi_Displaced wrt HLT_Dimuon0er16_Jpsi_NoOS_NoVertexing trig eff for Jpsi - "
jpsik4 = "HLT_DoubleMu4_3_Jpsi_Displaced wrt HLT_DoubleMu0 trig eff for Jpsi - "
jpsik5 = "L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4 trig eff for Jpsi - "
jpsik6 = "L1_DoubleMu0er1p6_dEta_Max1p8_OS trig eff for Jpsi - "
bsmm   = "HLT_DoubleMu4_3_Bs wrt HLT_Dimuon0_LowMass_L1_0er1p5 trig eff for Bs - "
bsmm2  = "HLT_DoubleMu4_3_Bs wrt HLT_DoubleMu0 trig eff for Bs - "
bsmm3  = "L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4 trig eff for Bs - "
bsmm4  = "L1_DoubleMu0er1p6_dEta_Max1p8_OS trig eff for Bs - "

table = [
    # ( "Run2016BF", (jpsik3 + "2016BF (channel: 0)", jpsik3 + "2016BF MC (channel: 0)"), ("", "") ),
    # ( "Run2016GH", (jpsik3 + "2016GH (channel: 0)", jpsik3 + "2016GH MC (channel: 0)"), ("", "") ),
    # ( "Run2016H",  (jpsik4 + "2016H (channel: 0)",  jpsik4 + "2016H MC (channel: 0)"),  (bsmm2 + "2016H (channel: 0)", bsmm2 + "2016H MC (channel: 0)") ),
    # ( "Run2017",   (jpsik1 + "2017 (channel: 0)",   jpsik1 + "2017 MC (channel: 0)"),   (bsmm + "2017 (channel: 0)",   bsmm + "2017 MC (channel: 0)") ),
    # ( "Run2018",   (jpsik2 + "2018 (channel: 0)",   jpsik2 + "2018 MC (channel: 0)"),   (bsmm + "2018 (channel: 0)",   bsmm + "2018 MC (channel: 0)") ),
    # ( "Run2016BF", (jpsik3 + "2016BF (channel: 1)", jpsik3 + "2016BF MC (channel: 1)"), ("", "") ),
    # ( "Run2016GH", (jpsik3 + "2016GH (channel: 1)", jpsik3 + "2016GH MC (channel: 1)"), ("", "") ),
    # ( "Run2016H",  (jpsik4 + "2016H (channel: 1)",  jpsik4 + "2016H MC (channel: 1)"),  (bsmm2 + "2016H (channel: 1)", bsmm2 + "2016H MC (channel: 1)") ),
    # ( "Run2017",   (jpsik1 + "2017 (channel: 1)",   jpsik1 + "2017 MC (channel: 1)"),   (bsmm + "2017 (channel: 1)",   bsmm + "2017 MC (channel: 1)") ),
    # ( "Run2018",   (jpsik2 + "2018 (channel: 1)",   jpsik2 + "2018 MC (channel: 1)"),   (bsmm + "2018 (channel: 1)",   bsmm + "2018 MC (channel: 1)") )
    ( "Run2016BF", (jpsik6 + "2016BF (channel: 0)", jpsik6 + "2016BF MC (channel: 0)"), (bsmm4 + "2016BF (channel: 0)", bsmm4 + "2016BF MC (channel: 0)") ), 
    ( "Run2016GH", (jpsik6 + "2016GH (channel: 0)", jpsik6 + "2016GH MC (channel: 0)"), (bsmm4 + "2016GH (channel: 0)", bsmm4 + "2016GH MC (channel: 0)") ), 
    ( "Run2017",   (jpsik5 + "2017 (channel: 0)",   jpsik5 + "2017 MC (channel: 0)"),   (bsmm3 + "2017 (channel: 0)",   bsmm3 + "2017 MC (channel: 0)") ), 
    ( "Run2018",   (jpsik5 + "2018 (channel: 0)",   jpsik5 + "2018 MC (channel: 0)"),   (bsmm3 + "2018 (channel: 0)",   bsmm3 + "2018 MC (channel: 0)") ), 
    ( "Run2016BF", (jpsik6 + "2016BF (channel: 1)", jpsik6 + "2016BF MC (channel: 1)"), (bsmm4 + "2016BF (channel: 1)", bsmm4 + "2016BF MC (channel: 1)") ), 
    ( "Run2016GH", (jpsik6 + "2016GH (channel: 1)", jpsik6 + "2016GH MC (channel: 1)"), (bsmm4 + "2016GH (channel: 1)", bsmm4 + "2016GH MC (channel: 1)") ), 
    ( "Run2017",   (jpsik5 + "2017 (channel: 1)",   jpsik5 + "2017 MC (channel: 1)"),   (bsmm3 + "2017 (channel: 1)",   bsmm3 + "2017 MC (channel: 1)") ), 
    ( "Run2018",   (jpsik5 + "2018 (channel: 1)",   jpsik5 + "2018 MC (channel: 1)"),   (bsmm3 + "2018 (channel: 1)",   bsmm3 + "2018 MC (channel: 1)") ), 
] 

def ratio(data, mc):
    r = data['eff'] / mc['eff']
    err = sqrt(pow(data['eff_err'] / data['eff'], 2) + pow(mc['eff_err'] / mc['eff'], 2))
    return (r, r * err)

def print_results(data, mc):
    if data in results:
        print "& $%5.2f \pm %4.2f$ " % (results[data]['eff']*100, results[data]['eff_err']*100),
    else:
        print "&                  ",

    if mc in results:
        print "& $%5.2f \pm %4.2f$ " % (results[mc  ]['eff']*100, results[mc  ]['eff_err']*100),
    else:
        print "&                  ",

    if data in results and mc in results:
        print "& $%5.3f \pm %5.3f$ " % ratio(results[data], results[mc  ]),
    else:
        print "&                  ",


for ( year, (jpsik_data, jpsik_mc), (bsmm_data, bsmm_mc) ) in table:
    print "%- 9s " % (year),
    print_results(jpsik_data, jpsik_mc)
    print_results(bsmm_data,  bsmm_mc)
    print r'\\'
