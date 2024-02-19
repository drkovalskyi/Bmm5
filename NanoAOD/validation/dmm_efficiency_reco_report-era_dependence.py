import ROOT
import glob

from base_efficiency_reco_report import EfficiencyReport

limit = 10000
path = "/eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/529"
samples = [
    # {
    #     'final_state':'kpi',
    #     'name':'\dzkpi 22',
    #     'files': glob.glob(f'{path}/DstarToD0Pi_D0ToKPi_KPiFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v2+MINIAODSIM/*.root')[:limit],
    #     'chain':None
    # },
    # {
    #     'final_state':'kpi',
    #     'name':'\dzkpi 22EE',
    #     'files': glob.glob(f'{path}/DstarToD0Pi_D0ToKPi_KPiFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2+MINIAODSIM/*.root')[:limit],
    #     'chain':None
    # },
    # {
    #     'final_state':'kpi',
    #     'name':'\dzkpi 23',
    #     'files': glob.glob(f'{path}/DstarToD0Pi_D0ToKPi_KPiFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v15-v1+MINIAODSIM/*.root')[:limit],
    #     'chain':None
    # },
    # {
    #     'final_state':'kpi',
    #     'name':'\dzkpi 23BPix',
    #     'files': glob.glob(f'{path}/DstarToD0Pi_D0ToKPi_KPiFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer23BPixMiniAODv4-130X_mcRun3_2023_realistic_postBPix_v5-v1+MINIAODSIM/*.root')[:limit],
    #     'chain':None
    # },
    {
        'final_state':'mm',
        'name':'\dzmm 22',
        'files': glob.glob(f'{path}/DstarToD0Pi_D0To2Mu_MuFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v1+MINIAODSIM/*.root')[:limit],
        'chain':None
    },
    {
        'final_state':'mm',
        'name':'\dzmm 22EE',
        'files': glob.glob(f'{path}/DstarToD0Pi_D0To2Mu_MuFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22EEMiniAODv3-124X_mcRun3_2022_realistic_postEE_v1-v1+MINIAODSIM/*.root')[:limit],
        'chain':None
    },
    {
        'final_state':'mm',
        'name':'\dzmm 23',
        'files': glob.glob(f'{path}/DstarToD0Pi_D0To2Mu_MuFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v14-v1+MINIAODSIM/*.root')[:limit],
        'chain':None
    },
    {
        'final_state':'mm',
        'name':'\dzmm 23BPix',
        'files': glob.glob(f'{path}/DstarToD0Pi_D0To2Mu_MuFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer23BPixMiniAODv4-130X_mcRun3_2023_realistic_postBPix_v2-v1+MINIAODSIM/*.root')[:limit],
        'chain':None
    },
    # {
    #     'final_state':'pipi',
    #     'name':'\dzpipi 22',
    #     'files': glob.glob(f'{path}/DstarToD0Pi_D0To2Pi_PiFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v1+MINIAODSIM/*.root')[:limit],
    #     'chain':None
    # },
    # {
    #     'final_state':'pipi',
    #     'name':'\dzpipi 22EE',
    #     'files': glob.glob(f'{path}/DstarToD0Pi_D0To2Pi_PiFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22EEMiniAODv3-124X_mcRun3_2022_realistic_postEE_v1-v1+MINIAODSIM/*.root')[:limit],
    #     'chain':None
    # },
    # {
    #     'final_state':'pipi',
    #     'name':'\dzpipi 23',
    #     'files': glob.glob(f'{path}/DstarToD0Pi_D0To2Pi_PiFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v14-v1+MINIAODSIM/*.root')[:limit],
    #     'chain':None
    # },
    # {
    #     'final_state':'pipi',
    #     'name':'\dzpipi 23BPix',
    #     'files': glob.glob(f'{path}/DstarToD0Pi_D0To2Pi_PiFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer23BPixMiniAODv4-130X_mcRun3_2023_realistic_postBPix_v2-v1+MINIAODSIM/*.root')[:limit],
    #     'chain':None
    # },
]

cuts = [
    {
        'cut':{
            'mm':'dstar_gen_pdgId!=0 && dstar_mm_index>=0 && mm_mu1_pt[dstar_mm_index]>4 && mm_mu2_pt[dstar_mm_index]>4',
            'pipi':(
                'dstar_hh_index>=0 && hh_had1_pdgId[dstar_hh_index]*hh_had2_pdgId[dstar_hh_index]==-211*211 && '
                'abs(dstar_gen_pdgId)==413 && hh_gen_had1_pdgId[dstar_hh_index]*hh_gen_had2_pdgId[dstar_hh_index]==-211*211 &&'
                'hh_had1_pt[dstar_hh_index]>4 && hh_had2_pt[dstar_hh_index]>4'
            ),
            'kpi':(
                'dstar_hh_index>=0 && hh_had1_pdgId[dstar_hh_index]*hh_had2_pdgId[dstar_hh_index]==-321*211 && '
                'abs(dstar_gen_pdgId)==413 && hh_gen_had1_pdgId[dstar_hh_index]*hh_gen_had2_pdgId[dstar_hh_index]==-321*211 &&'
                'hh_had1_pt[dstar_hh_index]>4 && hh_had2_pt[dstar_hh_index]>4'
            ),
            'fakes':'dstar_gen_pdgId!=0 && dstar_mm_index>=0 && abs(mm_gen_cpdgId[dstar_mm_index])==421 && mm_mu1_pt[dstar_mm_index]>4 && mm_mu2_pt[dstar_mm_index]>4',
        },
        'name':'MC matched baseline',
    },
    {
        'cut':{
            'mm':'HLT_DoubleMu4_3_LowMass',
            'fakes':'HLT_DoubleMu4_3_LowMass',
        },
        'name':r'\verb|HLT_DoubleMu4_3_LowMass|',
    },
    {
        'cut':{
            'mm':'dstar_dm_pv>0.140 && dstar_dm_pv<0.155',
            'pipi':'dstar_dm_pv>0.140 && dstar_dm_pv<0.155',
            'kpi':'dstar_dm_pv>0.140 && dstar_dm_pv<0.155',
            'fakes':'dstar_dm_pv>0.140 && dstar_dm_pv<0.155',
        },
        'name':'$\dm$ in [0.140, 0.155]',
    },
    {
        'cut':{
            'mm':'mm_kin_mass[dstar_mm_index]>1.81 && mm_kin_mass[dstar_mm_index]<1.94',
            'pipi':'hh_kin_mass[dstar_hh_index]>1.81 && hh_kin_mass[dstar_hh_index]<1.94',
            'kpi':'hh_kin_mass[dstar_hh_index]>1.81 && hh_kin_mass[dstar_hh_index]<1.94',
            'fakes':'mm_kin_mass[dstar_mm_index]>1.81 && mm_kin_mass[dstar_mm_index]<1.94',
        },
        'name':'$\PDz$ mass in [1.81, 1.94]',
    },
    {
        'cut':{
            'mm':'mm_kin_vtx_prob[dstar_mm_index]>0.01',
            'pipi':'hh_kin_vtx_prob[dstar_hh_index]>0.01',
            'kpi':'hh_kin_vtx_prob[dstar_hh_index]>0.01',
            'fakes':'mm_kin_vtx_prob[dstar_mm_index]>0.01',
        },
        'name':'$\PDz$ vertex probability $> 0.01$',
    },
    {
        'cut':{
            'mm':'dstar_pv_with_pion_prob>0.1',
            'pipi':'dstar_pv_with_pion_prob>0.1',
            'kpi':'dstar_pv_with_pion_prob>0.1',
            'fakes':'dstar_pv_with_pion_prob>0.1',
        },
        'name':'$\PDstpm$ vertex probability $> 0.1$',
    },
    # {
    #     'cut':{
    #         'mm':'Muon_softMva[mm_mu1_index[dstar_mm_index]] > 0.45 && Muon_softMva[mm_mu2_index[dstar_mm_index]] > 0.45',
    #         'fakes':'Muon_softMva[mm_mu1_index[dstar_mm_index]] > 0.45 && Muon_softMva[mm_mu2_index[dstar_mm_index]] > 0.45',
    #     },
    #     'name':'Muon identification',
    # },
    {
        'cut':{
            'mm':'mm_kin_sl3d[dstar_mm_index]>3',
            'pipi':'hh_kin_sl3d[dstar_hh_index]>3',
            'kpi':'hh_kin_sl3d[dstar_hh_index]>3',
            'fakes':'mm_kin_sl3d[dstar_mm_index]>3',
        },
        'name':r'$\fls > 3$',
    },
    {
        'cut':{
            'mm':'mm_kin_alpha[dstar_mm_index]<0.1',
            'pipi':'hh_kin_alpha[dstar_hh_index]<0.1',
            'kpi':'hh_kin_alpha[dstar_hh_index]<0.1',
            'fakes':'mm_kin_alpha[dstar_mm_index]<0.1',
        },
        'name':'$\pa < 0.1$',
    },
    
]

report = EfficiencyReport(samples, cuts)
report.make_report("gen", r"%8.5f")

