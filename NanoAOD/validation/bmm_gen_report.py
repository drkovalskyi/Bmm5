import ROOT, copy

samples = [
    {
        'type':'mm',
        'name':'\bsmm',
        'files':[
            "/eos/cms/store/group/phys_muon/dmytro/tmp/RunIISummer20UL18MiniAOD_BsToMuM_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen.root"
        ],
    },
    {
        'type':'mmk',
        'name':'\bjpsik',
        'files':[
            "/eos/cms/store/group/phys_muon/dmytro/tmp/RunIISummer20UL18MiniAOD_BuToJpsiK_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen.root"
        ],
    },
    {
        'type':'mmkk',
        'name':'blah',
        'files':[
            "/eos/cms/store/group/phys_muon/dmytro/tmp/RunIISummer20UL18MiniAOD_BsToJPsiPhi_JPsiToMuMu_PhiToKK_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen.root"
        ]
    }
]

cuts = [
    {
        'cut':'genbmm_mu1_pt>3.5&&genbmm_mu2_pt>3.5',
        'name':'$\ptmuone>3.5$, $\ptmutwo>3.5$',
        'types':['mm','mmk','mmkk']
    },
    {
        'cut':'abs(genbmm_mu1_eta)<2.5&&abs(genbmm_mu2_eta)<2.5',
        'name':'$|\etamuone|<2.5$, $|\etamutwo|<2.5|$',
        'types':['mm','mmk','mmkk']
    },
    {
        'cut':'genbmm_dau3_pt>0.4&&abs(genbmm_dau3_eta)<2.5',
        'name':'kaon $\pt>0.4$, $|\eta|<2.5$',
        'types':['mmk','mmkk']
    },
    {
        'cut':'genbmm_dau4_pt>0.4&&abs(genbmm_dau4_eta)<2.5',
        'name':'two kaons $\pt>0.4$, $|\eta|<2.5$',
        'types':['mmkk']
    },
    {
        'cut':'genbmm_mu1_pt>4.0&&genbmm_mu2_pt>4.0',
        'name':'$\ptmuone>4.0$, $\ptmutwo>4.0$',
        'types':['mm','mmk','mmkk']
    },
    {
        'cut':'abs(genbmm_mu1_eta)<1.4&&abs(genbmm_mu2_eta)<1.4',
        'name':'$|\etamuone|<1.4$, $|\etamutwo|<1.4$',
        'types':['mm','mmk','mmkk']
    },
    {
        'cut':'genbmm_dau3_pt>1.0',
        'name':'kaon $\pt>1.0$',
        'types':['mmk','mmkk']
    },
    {
        'cut':'genbmm_dau4_pt>1.0',
        'name':'two kaons $\pt>1.0$',
        'types':['mmkk']
    },
    {
        'cut':'genbmm_mu1_good && genbmm_mu2_good',
        'name':'Two good reco muons',
        'types':['mm','mmk','mmkk']
    },
    {
        'cut':'genbmm_dau3_reco_pt>0',
        'name':'Reco kaon',
        'types':['mmk','mmkk']
    },
    {
        'cut':'genbmm_dau4_reco_pt>0',
        'name':'Twi reco kaons',
        'types':['mmkk']
    },
]

chains = []
totals = []

for sample in samples:
    chain = ROOT.TChain("Events")
    for f in sample['files']:
        chain.Add(f)
    chains.append(chain)
    totals.append(chain.GetEntries())

effective_cuts = [""] * len(chains)
n_last = copy.deepcopy(totals)

for entry in cuts:
    print("%40s " % entry['name'], end=' ')
    for i,chain in enumerate(chains):
        if samples[i]['type'] in entry['types']:
            if effective_cuts[i] != "":  effective_cuts[i] += "&&"
            effective_cuts[i] += entry['cut']
            n = chain.GetEntries(effective_cuts[i])
            print("& %6.2f & %6.2f " % (100.0 * n / totals[i], 100.0 * n / n_last[i]), end=' ')
            n_last[i] = n
        else:
            print("& & ", end=' ')
        
    print("\\\\")
    
