import ROOT

samples = [
    {
        'final_state':'mm',
        'name':'\bsmm',
        'files':[
            "/eos/cms/store/group/phys_muon/dmytro/tmp/RunIISummer20UL18MiniAOD_BsToMuM_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen.root"
        ],
        'chain':None
    },
    {
        'final_state':'mmk',
        'name':'\bjpsik',
        'files':[
            "/eos/cms/store/group/phys_muon/dmytro/tmp/RunIISummer20UL18MiniAOD_BuToJpsiK_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen.root"
        ],
        'chain':None
    },
    {
        'final_state':'mmkk',
        'name':'blah',
        'files':[
            "/eos/cms/store/group/phys_muon/dmytro/tmp/RunIISummer20UL18MiniAOD_BsToJPsiPhi_JPsiToMuMu_PhiToKK_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen.root"
        ],
        'chain':None
    }
]

cuts = [
    {
        'cut':{
            'mm':'mm_mu1_index>=0 && mm_mu2_index>=0 && mm_gen_pdgId!=0 && mm_kin_mu1pt>4 && mm_kin_mu2pt>4 && abs(mm_kin_mu1eta)<1.4 && abs(mm_kin_mu2eta)<1.4',
            'mmk':'mm_mu1_index[bkmm_mm_index]>=0 && mm_mu2_index[bkmm_mm_index]>=0 && bkmm_gen_pdgId!=0'
            + ' && Muon_pt[mm_mu1_index[bkmm_mm_index]]>4 && Muon_pt[mm_mu2_index[bkmm_mm_index]]>4'
            + ' && abs(Muon_eta[mm_mu2_index[bkmm_mm_index]])<1.4 && abs(Muon_eta[mm_mu2_index[bkmm_mm_index]])<1.4',
            'mmkk':'mm_mu1_index[bkkmm_mm_index]>=0 && mm_mu2_index[bkkmm_mm_index]>=0 && bkkmm_gen_pdgId!=0'
            + ' && Muon_pt[mm_mu1_index[bkkmm_mm_index]]>4 && Muon_pt[mm_mu2_index[bkkmm_mm_index]]>4'
            + ' && abs(Muon_eta[mm_mu2_index[bkkmm_mm_index]])<1.4 && abs(Muon_eta[mm_mu2_index[bkkmm_mm_index]])<1.4',
        },
        'name':'MC matched baseline',
    },
    # {
    #     'cut':{
    #         'mm':'mm_kin_mu1pt>4 && mm_kin_mu2pt>4',
    #         'mmk':'Muon_pt[mm_mu1_index[bkmm_mm_index]]>4 && Muon_pt[mm_mu2_index[bkmm_mm_index]]>4',
    #         'mmkk':'Muon_pt[mm_mu1_index[bkkmm_mm_index]]>4 && Muon_pt[mm_mu2_index[bkkmm_mm_index]]>4',
    #     },
    #     'name':'$\ptmuone>4.0$,$\ptmutwo>4.0$',
    # },
    # {
    #     'cut':{
    #         'mm':'abs(mm_kin_mu1eta)<1.4 && abs(mm_kin_mu2eta)<1.4',
    #         'mmk':'abs(Muon_eta[mm_mu2_index[bkmm_mm_index]])<1.4 && abs(Muon_eta[mm_mu2_index[bkmm_mm_index]])<1.4',
    #         'mmkk':'abs(Muon_eta[mm_mu2_index[bkkmm_mm_index]])<1.4 && abs(Muon_eta[mm_mu2_index[bkkmm_mm_index]])<1.4',
    #     },
    #     'name':'$|\etamuone|<1.4$,$|\etatwo|<1.4$',
    # },
    {
        'cut':{
            'mm':'mm_kin_pt>5',
            'mmk':'bkmm_jpsimc_pt>5',
            'mmkk':'bkkmm_jpsikk_pt>5',
        },
        'name':'$\ptmm>5.0$',
    },
    {
        'cut':{
            'mmk':'bkmm_jpsimc_pt>7',
            'mmkk':'bkkmm_jpsikk_pt>7',
        },
        'name':'$\ptmm>7.0$',
    },
    {
        'cut':{
            'mmkk':'bkkmm_jpsikk_kk_mass>0.99&&bkkmm_jpsikk_kk_mass<1.05',
        },
        'name':'$\phi$ mass in [0.99,1.05]',
    },
    {
        'cut':{
            'mm':'abs(mm_kin_mass-5.4)<0.5',
            'mmk':'abs(bkmm_jpsimc_mass-5.4)<0.5',
            'mmkk':'abs(bkkmm_jpsikk_mass-5.4)<0.5',
        },
        'name':'$B$ mass in [4.9,5.9]',
    },
    {
        'cut':{
            'mm':'mm_kin_sl3d>4',
            'mmk':'bkmm_jpsimc_sl3d>4',
            'mmkk':'bkkmm_jpsikk_sl3d>4',
        },
        'name':'$\\fls>4$',
    },
    {
        'cut':{
            'mm':'mm_kin_sl3d>6',
        },
        'name':'$\\fls>6$',
    },
    {
        'cut':{
            'mm':'mm_kin_vtx_prob>0.025',
            'mmk':'bkmm_jpsimc_vtx_prob>0.025',
            'mmkk':'bkkmm_jpsikk_vtx_prob>0.025',
        },
        'name':'vertex probability $> 0.025$'
    },
    {
        'cut':{
            'mmk':'bkmm_jpsimc_vtx_prob>0.1',
            'mmkk':'bkkmm_jpsikk_vtx_prob>0.1',
        },
        'name':'vertex probability $> 0.1$'
    },
    {
        'cut':{
            'mm':'Muon_softMva[mm_mu1_index]>0.45&&Muon_softMva[mm_mu2_index]>0.45',
            'mmk':'Muon_softMva[mm_mu1_index[bkmm_mm_index]]>0.45&&Muon_softMva[mm_mu2_index[bkmm_mm_index]]>0.45',
            'mmkk':'Muon_softMva[mm_mu1_index[bkkmm_mm_index]]>0.45&&Muon_softMva[mm_mu2_index[bkkmm_mm_index]]>0.45',
        },
        'name':'soft MVA ID $>0.45$'
    },
    {
        'cut':{
            'mm':'HLT_DoubleMu4_3_Bs',
            'mmk':'HLT_DoubleMu4_3_Jpsi',
            'mmkk':'HLT_DoubleMu4_3_Jpsi',
        },
        'name':'trigger'
    },
    {
        'cut':{
            'mm':'mm_mva>0.90',
        },
        'name':'MVA > 0.90'
    },
    {
        'cut':{
            'mm':'mm_mva>0.99',
        },
        'name':'MVA > 0.99'
    },
    
]


total_counts = []
current_counts = []
final_counts = []

current_cuts = [""] * len(samples)

# function to get N and N-1 cuts
def get_final_cut(final_state, cut_to_exclude=None):
    cut = ""
    for entry in cuts:
        if final_state in entry['cut']:
            if cut_to_exclude and cut_to_exclude == entry['cut'][final_state]:
                continue
            if cut != "":  cut += "&&"
            cut += entry['cut'][final_state]
    return "mm_mu1_index>=0 && mm_mu2_index>=0 && " + cut

for sample in samples:
    chain = ROOT.TChain("Events")
    for f in sample['files']:
        chain.Add(f)
    sample['chain'] = chain

    n = chain.GetEntries()
    total_counts.append(n)
    current_counts.append(n)
    
    cut = get_final_cut(sample['final_state'])
    final_counts.append(chain.GetEntries(cut))


for entry in cuts:
    print "%35s " % entry['name'],
    for i,sample in enumerate(samples):
        if sample['final_state'] in entry['cut']:
            baseline_cut = True
            if current_cuts[i] != "":
                current_cuts[i] += "&&"
                baseline_cut = False
            current_cuts[i] += entry['cut'][sample['final_state']]
            n1 = sample['chain'].GetEntries(current_cuts[i])
            print "& %6.2f & %6.2f " % (100.0 * n1 / total_counts[i], 100.0 * n1 / current_counts[i]),
            current_counts[i] = n1
            if not baseline_cut:
                n2 = sample['chain'].GetEntries(get_final_cut(sample['final_state'], entry['cut'][sample['final_state']]))
                print "& %6.2f " % (100.0 * final_counts[i] / n2),
            else:
                print "&        ",
        else:
            print "&        &         &        ",
        
    print "\\\\"
    
