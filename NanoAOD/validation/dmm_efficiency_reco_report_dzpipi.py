import ROOT

from base_efficiency_reco_report import EfficiencyReport

path = "/eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/523"
path2 = "/eos/cms/store/group/phys_bphys/bmm/bmm6/PostProcessing/Skims/523"
samples = [
    {
        'final_state':'c_pipi',
        'name':'\dzpipi (c-quark)',
        'files':[
            "/eos/cms/store/group/phys_muon/dmytro/tmp/DstarToD0Pi_D0To2Pi_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen.root",
        ],
        'chain':None
    },
    {
        'final_state':'pipi',
        'name':'\dzpipi (minbias)',
        'files':[
            path2 + "/dzpipi/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8+Run3Summer22MiniAODv3-Pilot_124X_mcRun3_2022_realistic_v12-v5+MINIAODSIM/09df1c6dc6b9faa5443c7e7c40002838.root",
            path2 + "/dzpipi/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8+Run3Summer22MiniAODv3-Pilot_124X_mcRun3_2022_realistic_v12-v5+MINIAODSIM/2dc95a5b20db06295ee5c8beab54f825.root",
            path2 + "/dzpipi/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8+Run3Summer22MiniAODv3-Pilot_124X_mcRun3_2022_realistic_v12-v5+MINIAODSIM/3aae8fc447ea454ab7e452e123f17097.root",
            path2 + "/dzpipi/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8+Run3Summer22MiniAODv3-Pilot_124X_mcRun3_2022_realistic_v12-v5+MINIAODSIM/407920bde2864ed70503fc806d78f315.root",
            path2 + "/dzpipi/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8+Run3Summer22MiniAODv3-Pilot_124X_mcRun3_2022_realistic_v12-v5+MINIAODSIM/4f01c0c8ad757d923abf4b1f2a651bfc.root",
            path2 + "/dzpipi/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8+Run3Summer22MiniAODv3-Pilot_124X_mcRun3_2022_realistic_v12-v5+MINIAODSIM/520f6e55abdf40a053db4375c391a07a.root",
            path2 + "/dzpipi/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8+Run3Summer22MiniAODv3-Pilot_124X_mcRun3_2022_realistic_v12-v5+MINIAODSIM/5424a8672ca71aa5d83cdde55e647eb3.root",
            path2 + "/dzpipi/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8+Run3Summer22MiniAODv3-Pilot_124X_mcRun3_2022_realistic_v12-v5+MINIAODSIM/77f9cdf0d987b75fb8905b06e3675dde.root",
            path2 + "/dzpipi/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8+Run3Summer22MiniAODv3-Pilot_124X_mcRun3_2022_realistic_v12-v5+MINIAODSIM/7c41ef33106b41e82fb7fd3fd7e42c34.root",
            path2 + "/dzpipi/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8+Run3Summer22MiniAODv3-Pilot_124X_mcRun3_2022_realistic_v12-v5+MINIAODSIM/81b13389a89d47b8b9b7b4f781f5a297.root",
            path2 + "/dzpipi/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8+Run3Summer22MiniAODv3-Pilot_124X_mcRun3_2022_realistic_v12-v5+MINIAODSIM/8571965c33860c30f1e8de80ba618717.root",
            path2 + "/dzpipi/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8+Run3Summer22MiniAODv3-Pilot_124X_mcRun3_2022_realistic_v12-v5+MINIAODSIM/b3e80525b9ed48b4a74dae4a5b9fb7d0.root",
            path2 + "/dzpipi/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8+Run3Summer22MiniAODv3-Pilot_124X_mcRun3_2022_realistic_v12-v5+MINIAODSIM/e29a7ca7544e15299f871023ec07c012.root",
            path2 + "/dzpipi/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8+Run3Summer22MiniAODv3-Pilot_124X_mcRun3_2022_realistic_v12-v5+MINIAODSIM/ec70bf5453d852365f14704b85360909.root",
            path2 + "/dzpipi/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8+Run3Summer22MiniAODv3-Pilot_124X_mcRun3_2022_realistic_v12-v5+MINIAODSIM/fcadfeb58b2e32bb91f7f867a8440d79.root",
            path2 + "/dzpipi/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8+Run3Summer22MiniAODv3-Pilot_124X_mcRun3_2022_realistic_v12-v5+MINIAODSIM/fd2bdb6ecab3fcf8353dd112a72b0719.root",
            path2 + "/dzpipi/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8+Run3Summer22MiniAODv3-Pilot_124X_mcRun3_2022_realistic_v12-v5+MINIAODSIM/ff65166e7b41914bd632b2b130f0a6fc.root",

        ],
        'chain':None
    },
    {
        'final_state':'b_pipi',
        'name':'\dzpipi (b-quark)',
        'files':[
            path2 + "/dzpipi/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8+Run3Summer22MiniAODv3-Pilot_124X_mcRun3_2022_realistic_v12-v5+MINIAODSIM/09df1c6dc6b9faa5443c7e7c40002838.root",
            path2 + "/dzpipi/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8+Run3Summer22MiniAODv3-Pilot_124X_mcRun3_2022_realistic_v12-v5+MINIAODSIM/2dc95a5b20db06295ee5c8beab54f825.root",
            path2 + "/dzpipi/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8+Run3Summer22MiniAODv3-Pilot_124X_mcRun3_2022_realistic_v12-v5+MINIAODSIM/3aae8fc447ea454ab7e452e123f17097.root",
            path2 + "/dzpipi/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8+Run3Summer22MiniAODv3-Pilot_124X_mcRun3_2022_realistic_v12-v5+MINIAODSIM/407920bde2864ed70503fc806d78f315.root",
            path2 + "/dzpipi/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8+Run3Summer22MiniAODv3-Pilot_124X_mcRun3_2022_realistic_v12-v5+MINIAODSIM/4f01c0c8ad757d923abf4b1f2a651bfc.root",
            path2 + "/dzpipi/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8+Run3Summer22MiniAODv3-Pilot_124X_mcRun3_2022_realistic_v12-v5+MINIAODSIM/520f6e55abdf40a053db4375c391a07a.root",
            path2 + "/dzpipi/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8+Run3Summer22MiniAODv3-Pilot_124X_mcRun3_2022_realistic_v12-v5+MINIAODSIM/5424a8672ca71aa5d83cdde55e647eb3.root",
            path2 + "/dzpipi/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8+Run3Summer22MiniAODv3-Pilot_124X_mcRun3_2022_realistic_v12-v5+MINIAODSIM/77f9cdf0d987b75fb8905b06e3675dde.root",
            path2 + "/dzpipi/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8+Run3Summer22MiniAODv3-Pilot_124X_mcRun3_2022_realistic_v12-v5+MINIAODSIM/7c41ef33106b41e82fb7fd3fd7e42c34.root",
            path2 + "/dzpipi/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8+Run3Summer22MiniAODv3-Pilot_124X_mcRun3_2022_realistic_v12-v5+MINIAODSIM/81b13389a89d47b8b9b7b4f781f5a297.root",
            path2 + "/dzpipi/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8+Run3Summer22MiniAODv3-Pilot_124X_mcRun3_2022_realistic_v12-v5+MINIAODSIM/8571965c33860c30f1e8de80ba618717.root",
            path2 + "/dzpipi/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8+Run3Summer22MiniAODv3-Pilot_124X_mcRun3_2022_realistic_v12-v5+MINIAODSIM/b3e80525b9ed48b4a74dae4a5b9fb7d0.root",
            path2 + "/dzpipi/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8+Run3Summer22MiniAODv3-Pilot_124X_mcRun3_2022_realistic_v12-v5+MINIAODSIM/e29a7ca7544e15299f871023ec07c012.root",
            path2 + "/dzpipi/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8+Run3Summer22MiniAODv3-Pilot_124X_mcRun3_2022_realistic_v12-v5+MINIAODSIM/ec70bf5453d852365f14704b85360909.root",
            path2 + "/dzpipi/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8+Run3Summer22MiniAODv3-Pilot_124X_mcRun3_2022_realistic_v12-v5+MINIAODSIM/fcadfeb58b2e32bb91f7f867a8440d79.root",
            path2 + "/dzpipi/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8+Run3Summer22MiniAODv3-Pilot_124X_mcRun3_2022_realistic_v12-v5+MINIAODSIM/fd2bdb6ecab3fcf8353dd112a72b0719.root",
            path2 + "/dzpipi/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8+Run3Summer22MiniAODv3-Pilot_124X_mcRun3_2022_realistic_v12-v5+MINIAODSIM/ff65166e7b41914bd632b2b130f0a6fc.root",

        ],
        'chain':None
    },
]

cuts = [
    {
        'cut':{
            'c_pipi':(
                'dstar_hh_index>=0 && hh_had1_pdgId*hh_had2_pdgId==-211*211 && '
                'abs(dstar_gen_pdgId)==413 && abs(hh_gen_had1_pdgId)==211 && abs(hh_gen_had2_pdgId)==211 && '
                'hh_had1_pt[dstar_hh_index]>4 && hh_had2_pt[dstar_hh_index]>4'
            ),
            'pipi':(
                'dstar_hh_index>=0 && hh_had1_pdgId*hh_had2_pdgId==-211*211 && '
                'abs(dstar_gen_pdgId)==413 && abs(hh_gen_had1_pdgId)==211 && abs(hh_gen_had2_pdgId)==211 && '
                'hh_had1_pt[dstar_hh_index]>4 && hh_had2_pt[dstar_hh_index]>4'
            ),
            'b_pipi':(
                'dstar_hh_index>=0 && hh_had1_pdgId*hh_had2_pdgId==-211*211 && '
                '(abs(dstar_gen_mpdgId)==511 || abs(dstar_gen_mpdgId)==521 || abs(dstar_gen_mpdgId)==531) && '
                'abs(dstar_gen_pdgId)==413 && abs(hh_gen_had1_pdgId)==211 && abs(hh_gen_had2_pdgId)==211 && '
                'hh_had1_pt[dstar_hh_index]>4 && hh_had2_pt[dstar_hh_index]>4'
            ),
        },
        'name':'MC matched baseline',
    },
    # {
    #     'cut':{
    #         'mm':'mm_mu1_pt[dstar_mm_index] > 4 && mm_mu2_pt[dstar_mm_index] > 4',
    #         'pipi':'hh_had1_pt[dstar_hh_index] > 4 && hh_had2_pt[dstar_hh_index] > 4',
    #         'fakes':'mm_mu1_pt[dstar_mm_index] > 4 && mm_mu2_pt[dstar_mm_index] > 4',
    #     },
    #     'name':r'$\pt>4$',
    # },
    {
        'cut':{
            'c_pipi':'dstar_dm_pv>0.140 && dstar_dm_pv<0.155',
            'b_pipi':'dstar_dm_pv>0.140 && dstar_dm_pv<0.155',
            'pipi':'dstar_dm_pv>0.140 && dstar_dm_pv<0.155',
        },
        'name':'$\dm$ in [0.140, 0.155]',
    },
    {
        'cut':{
            'c_pipi':'hh_kin_mass[dstar_hh_index]>1.81 && hh_kin_mass[dstar_hh_index]<1.94',
            'b_pipi':'hh_kin_mass[dstar_hh_index]>1.81 && hh_kin_mass[dstar_hh_index]<1.94',
            'pipi':'hh_kin_mass[dstar_hh_index]>1.81 && hh_kin_mass[dstar_hh_index]<1.94',
        },
        'name':'$\PDz$ mass in [1.81, 1.94]',
    },
    {
        'cut':{
            'c_pipi':'hh_kin_vtx_prob[dstar_hh_index]>0.01',
            'b_pipi':'hh_kin_vtx_prob[dstar_hh_index]>0.01',
            'pipi':'hh_kin_vtx_prob[dstar_hh_index]>0.01',
        },
        'name':'$\PDz$ vertex probability $> 0.01$',
    },
    {
        'cut':{
            'c_pipi':'dstar_pv_with_pion_prob>0.1',
            'b_pipi':'dstar_pv_with_pion_prob>0.1',
            'pipi':'dstar_pv_with_pion_prob>0.1',
        },
        'name':'$\PDstpm$ vertex probability $> 0.1$',
    },
    {
        'cut':{
            'c_pipi':'hh_kin_sl3d[dstar_hh_index]>3',
            'b_pipi':'hh_kin_sl3d[dstar_hh_index]>3',
            'pipi':'hh_kin_sl3d[dstar_hh_index]>3',
        },
        'name':r'$\fls > 3$',
    },
    {
        'cut':{
            'c_pipi':'hh_kin_alpha[dstar_hh_index]<0.1',
            'b_pipi':'hh_kin_alpha[dstar_hh_index]<0.1',
            'pipi':'hh_kin_alpha[dstar_hh_index]<0.1',
        },
        'name':'$\pa < 0.1$',
    },
    
]

report = EfficiencyReport(samples, cuts)
report.make_report("first_cut")
