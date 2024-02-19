import ROOT
import glob

class EfficiencyReport:
    def __init__(self, samples, cuts):
        self.samples = samples
        self.cuts = cuts
        self.events = dict()
        self.lumis = dict()


    def get_first_cut(self, final_state):
        """ Get the first cut """

        if final_state in self.cuts[0]['cut']:
            return self.cuts[0]['cut'][final_state]
        else:
            return ""


    def get_complete_selection(self, final_state, cut_to_exclude=None):
        """ Get N and N-1 cuts """

        cut = ""
        for entry in self.cuts:
            if final_state in entry['cut']:
                if cut_to_exclude and cut_to_exclude == entry['cut'][final_state]:
                    continue
                if cut != "":  cut += "&&"
                cut += entry['cut'][final_state]
        return cut


    # def get_gen_filter_info(self, chain):
    #     n_gen = 0
    #     n_passed = 0
    def get_events(self, sample):
        if sample['name'] not in self.events:
            chain = ROOT.TChain("Events")
            for f in sample['files']:
                chain.Add(f)
            self.events[sample['name']] = chain

        return self.events[sample['name']]


    def get_lumi(self, sample):
        """ Get gen filter information """
        if sample['name'] not in self.lumis:
            if len(sample['files'])==0:
                raise Exception("No files are given for sample: %s" % sample['name'])

            f = ROOT.TFile.Open(sample['files'][0])

            n_gen_before_filter = 0
            n_gen_after_filter = 0

            if f.Get("LuminosityBlocks"):
                f.Close()

                lumis = ROOT.TChain("LuminosityBlocks")
                for f in sample['files']:
                    lumis.Add(f)

                for lumi in lumis:
                    n_gen_before_filter += lumi.GenFilter_numEventsTotal
                    n_gen_after_filter  += lumi.GenFilter_numEventsPassed

                self.lumis[sample['name']] = [n_gen_before_filter, n_gen_after_filter]

            elif f.Get("info"):
                f.Close()

                lumis = ROOT.TChain("info")
                for f in sample['files']:
                    lumis.Add(f)

                for lumi in lumis:
                    n_gen_before_filter += lumi.n_gen_all
                    n_gen_after_filter  += lumi.n_gen_passed

            self.lumis[sample['name']] = [n_gen_before_filter, n_gen_after_filter]
                
        return self.lumis[sample['name']]
    

    def make_report(self, baseline="sample", format="%6.2f"):

        formatted_empty_string = " " * len(format % 0.0)
        default_empty_string = " " * len("%5.1f" % 0.0)
        event_counts = []
        baseline_counts = []
        gen_passed_counts = []
        current_counts = []
        final_counts = []

        current_cuts = [""] * len(self.samples)
        
        for sample in self.samples:
            events = self.get_events(sample)
            lumi = self.get_lumi(sample)

            n_events_in_sample = events.GetEntries()
            event_counts.append(n_events_in_sample)

            cut = self.get_complete_selection(sample['final_state'])
            final_counts.append(events.GetEntries(cut))

            if baseline == "sample":
                baseline_counts.append(n_events_in_sample)
            elif baseline == "gen":
                baseline_counts.append(lumi[0])
                gen_passed_counts.append(lumi[1])
            elif baseline == "first_cut":
                cut = self.get_first_cut(sample['final_state'])
                baseline_counts.append(events.GetEntries(cut))
            else:
                raise Exception("Unsupported report baseline: %s" % baseline)
            
            current_counts.append(baseline_counts[-1])

        text_length = 30

        for cut in self.cuts:
            if len(cut['name']) >= text_length:
                text_length = len(cut['name']) + 1

        prefix = "%%%us " % text_length

        if baseline == "gen":
            print(prefix % "Generator filter", end=' ')
            for i,sample in enumerate(self.samples):
                scale = 100.0
                if 'scale' in sample:
                    scale = sample['scale']
                gen_eff = scale * gen_passed_counts[i] / baseline_counts[i]
                print(("& " + format) % gen_eff, end='')
                print(f" & {default_empty_string} & {default_empty_string} ", end=' ')
            print("\\\\")

        first_line = True
        for icut,entry in enumerate(self.cuts):
            if baseline != "first_cut" or icut != 0:
                print(prefix % entry['name'], end=' ')
            for i,sample in enumerate(self.samples):
                if baseline == "first_cut" and icut == 0:
                    current_cuts[i] = entry['cut'][sample['final_state']]
                    continue
                chain = self.get_events(sample)
                if sample['final_state'] in entry['cut']:
                    scale = 100.0
                    if 'scale' in sample:
                        scale = sample['scale']
                    if current_cuts[i] != "":
                        current_cuts[i] += "&&"
                    current_cuts[i] += entry['cut'][sample['final_state']]

                    #print "cut ",current_cuts[i]
                    n1 = chain.GetEntries(current_cuts[i])

                    print(("& " + format) % (scale * n1 / baseline_counts[i]), end=' ')

                    if not first_line:
                        n2 = chain.GetEntries(self.get_complete_selection(sample['final_state'], entry['cut'][sample['final_state']]))
                        print(("& %5.1f & %5.1f ") % (100.0 * n1 / current_counts[i], 100.0 * final_counts[i] / n2), end=' ')
                    else:
                        print(f"& {default_empty_string} & {default_empty_string} ", end=' ')

                    current_counts[i] = n1
                else:
                    print(f"& {formatted_empty_string} & {default_empty_string} & {default_empty_string} ", end=' ')

            if baseline != "first_cut" or icut != 0:
                print("\\\\")
                first_line = False

if __name__ == "__main__":
    path = "/eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/526"

    samples = [
        # {
        #     'final_state':'mm',
        #     'name':'\dzmm',
        #     'scale':1000.,
        #     'files': [
        #         # f for f in glob.glob(path + "/DstarToD0Pi_D0To2Mu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/*.root")
        #         f for f in glob.glob(path + "/DstarToD0Pi_D0To2Mu_MuFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v1+MINIAODSIM/*.root")
        #     ],
        # },
        {
            'final_state':'pipi',
            'name':'\dzpipi',
            'scale':1000.,
            'files':[
                # "/eos/cms/store/group/phys_bphys/bmm/dmytro/tmp/DstarToD0Pi_D0To2Pi_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen.root",
                # "/afs/cern.ch/work/d/dmytro/projects/Run3-Bmm-NanoAODv12/src/DstarToD0Pi_D0To2Pi_b-quark_NANOAODSIM_Bmm.root"
                f for f in glob.glob(path + "/DstarToD0Pi_D0To2Pi_PiFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v1+MINIAODSIM/*.root")
            ],
        },
        # {
        #     'final_state':'mm',
        #     'name':'\dzpipimm',
        #     'scale':1e9/50/50,
        #     'files':[
        #         f for f in glob.glob(path + "/DstarToD0Pi_D0To2Pi_PiToMu_PiFilter_PiLifetime0p1_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v2+MINIAODSIM/*.root")
        #         f for f in glob.glob(path + "/DstarToD0Pi_D0To2Pi_PiToMu_PiFilter_PiLifetime0p02_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v3+MINIAODSIM/*.root")
        #     ],
        # },
    ]

    cuts = [
        {
            'cut':{
                'mm':'dstar_mm_index>=0 && abs(dstar_gen_pdgId)==413 && mm_mu1_pt[dstar_mm_index]>4 && mm_mu2_pt[dstar_mm_index]>4 && Muon_isGlobal[mm_mu1_index[dstar_mm_index]] && Muon_isGlobal[mm_mu2_index[dstar_mm_index]]',
                'pipi':'dstar_hh_index>=0&&abs(dstar_gen_pdgId)==413 && abs(hh_gen_had1_pdgId[dstar_hh_index])==211&&abs(hh_gen_had2_pdgId[dstar_hh_index])==211 && hh_had1_pt[dstar_hh_index]>4 && hh_had2_pt[dstar_hh_index]>4 && abs(hh_gen_pdgId[dstar_hh_index])==421 && abs(hh_had1_pdgId[dstar_hh_index])==211 && abs(hh_had2_pdgId[dstar_hh_index])==211',
            },
            'name':'MC matched preselection',
        },
        {
            'cut':{
                'mm':'HLT_DoubleMu4_3_LowMass',
            },
            'name':r'\verb|HLT_DoubleMu4_3_LowMass|',
        },
        {
            'cut':{
                'mm':'dstar_dm_pv>0.140 && dstar_dm_pv<0.155',
                'pipi':'dstar_dm_pv>0.140 && dstar_dm_pv<0.155',
            },
            'name':'$\dm$ in [0.140, 0.155]',
        },
        {
            'cut':{
                'mm':'mm_kin_mass[dstar_mm_index]>1.81 && mm_kin_mass[dstar_mm_index]<1.94',
                'pipi':'hh_kin_mass[dstar_hh_index]>1.81 && hh_kin_mass[dstar_hh_index]<1.94',
            },
            'name':'$\PDz$ mass in [1.81, 1.94]',
        },
        {
            'cut':{
                'mm':'mm_kin_vtx_prob[dstar_mm_index]>0.01',
                'pipi':'hh_kin_vtx_prob[dstar_hh_index]>0.01',
            },
            'name':'$\PDz$ vertex probability $> 0.01$',
        },
        {
            'cut':{
                'mm':'dstar_pv_with_pion_prob>0.1',
                'pipi':'dstar_pv_with_pion_prob>0.1',
            },
            'name':'$\PDstpm$ vertex probability $> 0.1$',
        },
        # {
        #     'cut':{
        #         'mm':'Muon_softMva[mm_mu1_index[dstar_mm_index]] > 0.45 && Muon_softMva[mm_mu2_index[dstar_mm_index]] > 0.45',
        #     },
        #     'name':'Muon identification',
        # },
        {
            'cut':{
                'mm':'mm_kin_sl3d[dstar_mm_index]>3',
                'pipi':'hh_kin_sl3d[dstar_hh_index]>3',
            },
            'name':r'$\fls > 3$',
        },
        {
            'cut':{
                'mm':'mm_kin_alpha[dstar_mm_index]<0.1',
                'pipi':'hh_kin_alpha[dstar_hh_index]<0.1',
            },
            'name':'$\pa < 0.1$',
        },
    ]

    report = EfficiencyReport(samples, cuts)

    # report.make_report()
    # print()
    report.make_report("gen", r"%7.4f")


