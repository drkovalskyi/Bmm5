import ROOT, copy

path1 = "/eos/cms/store/group/phys_bphys/bmm/bmm5/PostProcessing/FlatNtuples/516/fit/"
path2 = "/eos/cms/store/group/phys_bphys/bmm/bmm5/PostProcessing/FlatNtuples/516/fit-bkmm/"

eras = ['RunIISummer20UL16MiniAODAPV', 'RunIISummer20UL16MiniAOD', 'RunIISummer20UL17MiniAOD', 'RunIISummer20UL18MiniAOD']
types = ['fs', 'fu']

samples = [
    {
        'type':'fs',
        'era':'RunIISummer20UL18MiniAOD',
        'filtered':False,
        'files':[
            path1 + "BsToMuMu_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v2+MINIAODSIM/*.root"
        ],
        'chain':None
    },
    {
        'type':'fs',
        'era':'RunIISummer20UL18MiniAOD',
        'filtered':True,
        'files':[
            path1 + "BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v1+MINIAODSIM/*.root"
        ],
        'chain':None
    },
    {
        'type':'fu',
        'era':'RunIISummer20UL18MiniAOD',
        'filtered':False,
        'files':[
            path2 + "BuToJpsiK_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v2+MINIAODSIM/*.root"
        ],
        'chain':None
    },
    {
        'type':'fu',
        'era':'RunIISummer20UL18MiniAOD',
        'filtered':True,
        'files':[
            path2 + "BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v2+MINIAODSIM/*.root"
        ],
        'chain':None
    },

    {
        'type':'fs',
        'era':'RunIISummer20UL17MiniAOD',
        'filtered':False,
        'files':[
            path1 + "BsToMuMu_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL17MiniAOD-106X_mc2017_realistic_v6-v2+MINIAODSIM/*.root"
        ],
        'chain':None
    },
    {
        'type':'fs',
        'era':'RunIISummer20UL17MiniAOD',
        'filtered':True,
        'files':[
            path1 + "BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL17MiniAOD-106X_mc2017_realistic_v6-v1+MINIAODSIM/*.root"
        ],
        'chain':None
    },
    {
        'type':'fu',
        'era':'RunIISummer20UL17MiniAOD',
        'filtered':False,
        'files':[
            path2 + "BuToJpsiK_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL17MiniAOD-106X_mc2017_realistic_v6-v2+MINIAODSIM/*.root",
        ],
        'chain':None
    },
    {
        'type':'fu',
        'era':'RunIISummer20UL17MiniAOD',
        'filtered':True,
        'files':[
            path2 + "BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL17MiniAOD-106X_mc2017_realistic_v6-v2+MINIAODSIM/*.root"
        ],
        'chain':None
    },

    {
        'type':'fs',
        'era':'RunIISummer20UL16MiniAOD',
        'filtered':False,
        'files':[
            path1 + "BsToMuMu_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL16MiniAOD-106X_mcRun2_asymptotic_v13-v2+MINIAODSIM/*.root"
        ],
        'chain':None
    },
    {
        'type':'fs',
        'era':'RunIISummer20UL16MiniAOD',
        'filtered':True,
        'files':[
            path1 + "BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL16MiniAOD-106X_mcRun2_asymptotic_v13-v1+MINIAODSIM/*.root"
        ],
        'chain':None
    },
    {
        'type':'fu',
        'era':'RunIISummer20UL16MiniAOD',
        'filtered':False,
        'files':[
            path2 + "BuToJpsiK_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL16MiniAOD-106X_mcRun2_asymptotic_v13-v2+MINIAODSIM/*.root",
        ],
        'chain':None
    },
    {
        'type':'fu',
        'era':'RunIISummer20UL16MiniAOD',
        'filtered':True,
        'files':[
            path2 + "BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL16MiniAOD-106X_mcRun2_asymptotic_v13-v1+MINIAODSIM/*.root"
        ],
        'chain':None
    },

    {
        'type':'fs',
        'era':'RunIISummer20UL16MiniAODAPV',
        'filtered':False,
        'files':[
            path1 + "BsToMuMu_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2+MINIAODSIM/*.root"
        ],
        'chain':None
    },
    {
        'type':'fs',
        'era':'RunIISummer20UL16MiniAODAPV',
        'filtered':True,
        'files':[
            path1 + "BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v1+MINIAODSIM/*.root"
        ],
        'chain':None
    },
    {
        'type':'fu',
        'era':'RunIISummer20UL16MiniAODAPV',
        'filtered':False,
        'files':[
            path2 + "BuToJpsiK_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2+MINIAODSIM/*.root",
        ],
        'chain':None
    },
    {
        'type':'fu',
        'era':'RunIISummer20UL16MiniAODAPV',
        'filtered':True,
        'files':[
            path2 + "BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v1+MINIAODSIM/*.root"
        ],
        'chain':None
    },
    
]


for sample in samples:
    chain = ROOT.TChain("info")
    for f in sample['files']:
        chain.Add(f)
    sample['chain'] = chain

for era in eras:
    print "%35s " % era,
    eff_s = None
    eff_u = None
    for type in types:
        eff_filtered = None
        eff_unfiltered = None
        for sample in samples:
            if sample['era'] != era: continue
            if sample['type'] != type: continue
            n_gen_passed = 0
            n_gen_all = 0
            for entry in sample['chain']:
                n_gen_passed += entry.n_gen_passed
                n_gen_all += entry.n_gen_all
            if sample['filtered']:
                eff_filtered = float(n_gen_passed)/n_gen_all
            else:
                eff_unfiltered = float(n_gen_passed)/n_gen_all
        print "& %6.3f & %6.3f " % (100.0 * eff_unfiltered, 100.0 * eff_filtered),
        if type == 'fs': eff_s = eff_unfiltered
        if type == 'fu': eff_u = eff_unfiltered
            
    print " & %0.4f\\\\" % (eff_s/eff_u)
    
