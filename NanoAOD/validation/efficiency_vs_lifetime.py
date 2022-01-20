import os
import ROOT
import tdrstyle

tdrstyle.setTDRStyle()
ROOT.gROOT.SetBatch(True)
c1 = ROOT.TCanvas("c1","c1", 1200, 1200)

baseline = "mm_gen_pdgId!=0"

selection = "mm_gen_pdgId!=0 && "\
    "mm_mu1_index>=0 && mm_mu2_index>=0 && "\
    "Muon_softMva[mm_mu1_index]>0.45 && "\
    "abs(mm_kin_mu1eta)<1.4 && "\
    "mm_kin_mu1pt>4 && "\
    "Muon_softMva[mm_mu2_index]>0.45 && "\
    "abs(mm_kin_mu2eta)<1.4 && "\
    "mm_kin_mu2pt>4 && "\
    "abs(mm_kin_mass-5.4)<0.5 && "\
    "mm_kin_sl3d>6 && "\
    "mm_kin_pt>5.0 && mm_kin_vtx_prob>0.025"

mva_cut_loose = "mm_mva>0.90"
mva_cut_tight = "mm_mva>0.99"

samples = [
    {
        'lifetime':'1.40',
        'files':'/eos/cms/store/group/phys_bphys/bmm/bmm5/NanoAOD/516/BsToMuMu_SoftQCDnonD_TuneCP5_BsLifetime1p40_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v1+MINIAODSIM/*.root'
    },
    {
        'lifetime':'1.60',
        'files':'/eos/cms/store/group/phys_bphys/bmm/bmm5/NanoAOD/516/BsToMuMu_SoftQCDnonD_TuneCP5_BsLifetime1p60_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v1+MINIAODSIM/*.root'
    },
]

effs = dict()

def process(sample):
    lf = sample['lifetime'] 
    print "Processing sample with lifetime: %s" % lf
    effs[lf] = []
    
    lumis = ROOT.TChain("LuminosityBlocks")
    lumis.Add(sample['files'])
    events = ROOT.TChain("Events")
    events.Add(sample['files'])

    n_gen = 0
    n_passed = 0
    for lumi in lumis:
        n_gen += lumi.GenFilter_numEventsTotal
        n_passed += lumi.GenFilter_numEventsPassed

    print "Gen efficiency: %0.3f" % (100. * n_passed / n_gen)

    n_baseline = events.GetEntries(baseline)
    print "Baseline efficiency wrt passed: %0.3f" % (100. * n_baseline / n_passed)
    effs[lf].append(float(n_baseline) / n_passed)
    
    n_preselection = events.GetEntries(baseline + "&&" + selection)
    print "Preselection efficiency wrt passed: %0.3f" % (100. * n_preselection / n_passed)
    effs[lf].append(float(n_preselection) / n_passed)

    n_full_loose = events.GetEntries(baseline + "&&" + selection + "&&" + mva_cut_loose)
    print "MVA>0.90 efficiency wrt passed: %0.3f" % (100. * n_full_loose / n_passed)
    effs[lf].append(float(n_full_loose) / n_passed)

    n_full_tight = events.GetEntries(baseline + "&&" + selection + "&&" + mva_cut_tight)
    print "MVA>0.99 efficiency wrt passed: %0.3f" % (100. * n_full_tight / n_passed)
    effs[lf].append(float(n_full_tight) / n_passed)
    

for sample in samples:
    process(sample)
    print

for i, val in enumerate(effs['1.40']):
    print "%0.1f%%" % (100.*val/effs['1.60'][i])
