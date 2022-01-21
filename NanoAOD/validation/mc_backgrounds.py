import os, json
import ROOT
import tdrstyle
from pprint import pprint
from collections import OrderedDict
from math import *

bbbar_cross_section = 0.5343 #mb McM

lumi = 140. #1/fb

recompute_results = False
result_file = "results/mc_backgrounds.json"
results = dict()
if os.path.exists(result_file):
    results = json.load(open(result_file))

ROOT.gROOT.SetBatch(True)
c1 = ROOT.TCanvas("c1","c1", 800, 800)

# mc_match = "mm_gen_pdg!=0"

preselection_noid = "abs(mm_kin_mu1eta)<1.4 && "\
    "mm_kin_mu1pt>4 && "\
    "abs(mm_kin_mu2eta)<1.4 && "\
    "mm_kin_mu2pt>4 && "\
    "abs(mm_kin_mass-5.4)<0.5 && "\
    "mm_kin_sl3d>6 && "\
    "mm_kin_pt>5.0 && mm_kin_vtx_prob>0.025"

preselection = preselection_noid + "&&"\
    "mm_mu1_index>=0 && mm_mu2_index>=0 && HLT_DoubleMu4_3_Bs"

prefit = preselection + "&&"\
    "Muon_softMva[mm_mu1_index]>0.45 && "\
    "Muon_softMva[mm_mu2_index]>0.45"

mva050 = prefit + "&& mm_mva>0.50"
mva090 = prefit + "&& mm_mva>0.90"
final  = prefit + "&& mm_mva>0.99"
noid_mva099 = preselection_noid + "&& mm_mva>0.99"

path = "/eos/cms/store/group/phys_bphys/bmm/bmm5/NanoAOD/517/"

def get_ngen(files):
    lumis = ROOT.TChain("LuminosityBlocks")
    for entry in files:
        lumis.Add(entry)
    n = 0
    for lumi in lumis:
        n += lumi.GenFilter_numEventsTotal
    return n

def get_info(files, selection):
    events = ROOT.TChain("Events")
    for entry in files:
        events.Add(entry)
    hist = ROOT.TH1F("hist","",100,4.9,5.9)
    events.Draw("mm_kin_mass>>hist", selection)
    hist.SetDirectory(0)
    return {
        'n':hist.GetEntries(),
        'mean':hist.GetMean(),
        'rms':hist.GetRMS()
    }

def get_relative_yield(name, selection):
    bf_ratio = float(samples[name]['bf']) / samples[reference]['bf']
    eff      = float(results[name][selection]['n']) / results[name]['ngen']
    eff_err  = sqrt(results[name][selection]['n']) / results[name]['ngen']
    eff_ref  = float(results[reference][selection]['n']) / results[reference]['ngen']

    return bf_ratio * eff / eff_ref, bf_ratio * eff_err / eff_ref

def get_absolute_yield(name, selection):
    c = bbbar_cross_section * lumi * 1e12
    eff      = float(results[name][selection]['n']) / results[name]['ngen']
    eff_err  = sqrt(results[name][selection]['n']) / results[name]['ngen']

    return c * eff * samples[name]['bf'], c * eff_err * samples[name]['bf']


        
samples = {
    'BsToMuMu': {
        'files': [path + '/BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v1+MINIAODSIM/1*.root'],
        'bf': 3.66E-9,
        'name':r'B_s\to\mu\mu',
    },
    'BdToMuMu': {
        'files': [path + '/BdToMuMu_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v2+MINIAODSIM/*.root'],
        'bf': 1.03E-10,
        'name':r'B_d\to\mu\mu',
    },
    'BTohh': {
        'files':[
            path + '/BTohh_hToMuNu_BsBdMixture_modHadLifetime_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL16MiniAOD-106X_mcRun2_asymptotic_v13-v3+MINIAODSIM/*.root',
            path + '/BTohh_hToMuNu_BsBdMixture_modHadLifetime_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v3+MINIAODSIM/*.root',
            path + '/BTohh_hToMuNu_BsBdMixture_modHadLifetime_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL17MiniAOD-106X_mc2017_realistic_v6-v3+MINIAODSIM/*.root',
            path + '/BTohh_hToMuNu_BsBdMixture_modHadLifetime_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v3+MINIAODSIM/*.root',

        ],
        'bf': 3e-5/100,
        'name':r'B_x \to hh',
    },
    # 'LambdaBToPK': {
    #     'files':[
    #         path + '/LambdaBToPK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v1+MINIAODSIM/*.root',
    #         path + '/LambdaBToPK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL16MiniAOD-106X_mcRun2_asymptotic_v13-v1+MINIAODSIM/*.root',
    #         path + '/LambdaBToPK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v1+MINIAODSIM/*.root',
    #         path + '/LambdaBToPK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL17MiniAOD-106X_mc2017_realistic_v6-v1+MINIAODSIM/*.root'
    #     ],
    #     'bf': 5.4e-6,
    #     'name':r'\Lambda_B\to p K',
    # },
    # 'LambdaBToPPi': {
    #     'files':[
    #         path + '/LambdaBToPPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL16MiniAOD-106X_mcRun2_asymptotic_v13-v1+MINIAODSIM/*.root',
    #         path + '/LambdaBToPPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v1+MINIAODSIM/*.root',
    #         path + '/LambdaBToPPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL17MiniAOD-106X_mc2017_realistic_v6-v1+MINIAODSIM/*.root',
    #         path + '/LambdaBToPPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v1+MINIAODSIM/*.root',
    #     ],
    #     'bf': 4.5e-6,
    #     'name':r'\Lambda_B\to p \pi',
    # },
    # 'LambdaBToPMuNu': {
    #     'files':[
    #         path + '/LambdaBToPMuNu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL16MiniAOD-106X_mcRun2_asymptotic_v13-v1+MINIAODSIM/*.root',
    #         path + '/LambdaBToPMuNu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v1+MINIAODSIM/*.root',
    #         path + '/LambdaBToPMuNu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL17MiniAOD-106X_mc2017_realistic_v6-v1+MINIAODSIM/*.root',
    #     ],
    #     'bf': 4.1e-4,
    #     'name':r'\Lambda_B\to p \mu\nu',
    # },
    'BdToKK': {
        'files':[
            path + '/BdToKK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL16MiniAOD-106X_mcRun2_asymptotic_v13-v1+MINIAODSIM/*.root',
            path + '/BdToKK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v1+MINIAODSIM/*.root',
            path + '/BdToKK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL17MiniAOD-106X_mc2017_realistic_v6-v1+MINIAODSIM/*.root',
            path + '/BdToKK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v1+MINIAODSIM/*.root',
        ],
        'bf': 7.8e-8,
        'name':r'B_d\to K K',
    },
    'BsToKK':{
        'files':[
            path + '/BsToKK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL16MiniAOD-106X_mcRun2_asymptotic_v13-v1+MINIAODSIM/*.root',
            path + '/BsToKK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v1+MINIAODSIM/*.root',
            path + '/BsToKK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL17MiniAOD-106X_mc2017_realistic_v6-v1+MINIAODSIM/*.root',
            path + '/BsToKK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v1+MINIAODSIM/*.root',
        ],
        'bf': 2.66e-5,
        'name':r'B_s\to K K',
    },
    'BdToPiPi':{
        'files':[
            path + '/BdToPiPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL16MiniAOD-106X_mcRun2_asymptotic_v13-v1+MINIAODSIM/*.root',
            path + '/BdToPiPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v1+MINIAODSIM/*.root',
            path + '/BdToPiPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL17MiniAOD-106X_mc2017_realistic_v6-v1+MINIAODSIM/*.root',
            path + '/BdToPiPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v1+MINIAODSIM/*.root',
        ],
        'bf': 5.12e-6,
        'name':r'B_d\to \pi\pi',
    },
    'BsToPiPi':{
        'files':[
            path + '/BsToPiPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL16MiniAOD-106X_mcRun2_asymptotic_v13-v1+MINIAODSIM/*.root',
            path + '/BsToPiPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v1+MINIAODSIM/*.root',
            path + '/BsToPiPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL17MiniAOD-106X_mc2017_realistic_v6-v1+MINIAODSIM/*.root',
            path + '/BsToPiPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v1+MINIAODSIM/*.root',
        ],
        'bf': 7e-7,
        'name':r'B_s\to \pi\pi',
    },
    'BdToKPi':{
        'files':[
            path + '/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL16MiniAOD-106X_mcRun2_asymptotic_v13-v1+MINIAODSIM/*.root',
            path + '/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v1+MINIAODSIM/*.root',
            path + '/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL17MiniAOD-106X_mc2017_realistic_v6-v1+MINIAODSIM/*.root',
            path + '/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v1+MINIAODSIM/*.root',
        ],
        'bf':1.96e-5,
        'name':r'B_d\to K\pi',
    },
    'BsToKPi':{
        'files':[
            path + '/BsToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL16MiniAOD-106X_mcRun2_asymptotic_v13-v1+MINIAODSIM/*.root',
            path + '/BsToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v1+MINIAODSIM/*.root',
            path + '/BsToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL17MiniAOD-106X_mc2017_realistic_v6-v1+MINIAODSIM/*.root',
            path + '/BsToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v1+MINIAODSIM/*.root',
        ],
        'bf':5.8e-6,
        'name':r'B_s\to K\pi',
    }
    
}

selections = OrderedDict([('loose no id', preselection_noid), ('mva099 no id',noid_mva099) ,('loose',preselection), ('medium',prefit), ('mva050',mva050), ('mva090',mva090),  ('tight',final) ])

for name,sample in samples.items():
    print name
    if name in results and not recompute_results:
        print "Results are already available. Skip"
        continue

    ngen = get_ngen(sample['files'])
    print "n_gen:", ngen
    results[name] = {'ngen': ngen }
    
    for selection_name, selection in selections.items():
        info = get_info(sample['files'], selection)
        print "%s selection:" % selection_name
        print "\tn:", info['n']
        print "\tmass mean: %5.2f GeV" % info['mean']
        print "\tmass RMS:  %5.2f GeV" % info['rms']

        results[name][selection_name] = info

    json.dump(results, open(result_file, "w"))

# pprint(results)

reference = 'BsToMuMu'

# print "Reference:", reference
# print "Expected relative yield with respect to the reference:"
# for name,sample in sorted(samples.items()):
#     if name == reference: continue
#     print name

#     for selection_name, selection in selections.items():
#         eff, eff_err = get_relative_yield(name, selection_name)
#         print "\t%-10s selection: %0.2f \pm %0.2f" % (selection_name, 100. * eff, 100. * eff_err)

# for name,sample in sorted(samples.items()):
#     if name == reference: continue
#     print "$%30s$" % sample['name'], 

#     for selection_name in ('loose', 'medium', 'tight'):
#         eff, eff_err = get_relative_yield(name, selection_name)
#         print "& $%6.2f \pm %5.2f \%%$ " % (100. * eff, 100. * eff_err),

#     print r'\\'

print
for name,sample in sorted(samples.items()):
    print "$%30s$" % sample['name'], 

    for selection_name in ('loose', 'medium', 'tight'):
        n, n_err = get_absolute_yield(name, selection_name)
        print "& $%6.1f \pm %5.1f \%%$ " % (n, n_err),

    print r'\\'

print "Effective lumi per sample:"

for name,sample in sorted(samples.items()):
    print "%-30s" % name,
    c = bbbar_cross_section * lumi * 1e12
    print "%10.2f/fb" % (float(results[name]['ngen'])/(bbbar_cross_section * 1e12 * samples[name]['bf']))
    
