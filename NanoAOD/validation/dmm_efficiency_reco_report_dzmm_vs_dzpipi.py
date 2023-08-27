import ROOT

from base_efficiency_reco_report import EfficiencyReport

path = "/eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/523"
path2 = "/eos/cms/store/group/phys_bphys/bmm/bmm6/PostProcessing/Skims/523"
samples = [
    {
        'final_state':'mm',
        'name':'\dzmm',
        'files':[
            path + "/DstarToD0Pi_D0To2Mu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/089fcee0-0260-470b-8f08-a458129f2c4a.root",
            path + "/DstarToD0Pi_D0To2Mu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/164a5846-abdf-44ae-ba68-00ac4669fdff.root",
            path + "/DstarToD0Pi_D0To2Mu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/1fbe986c-bc25-4a90-a0ea-03d2d4f76001.root",
            path + "/DstarToD0Pi_D0To2Mu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/27fb3f48-a6ea-4674-9fdd-af87426e0c0d.root",
            path + "/DstarToD0Pi_D0To2Mu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/2f2337ac-420c-477c-a781-1065f0c218a8.root",
            path + "/DstarToD0Pi_D0To2Mu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/2ff197b6-51fc-404f-8109-1db55acddc01.root",
            path + "/DstarToD0Pi_D0To2Mu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/3236c958-eb94-4338-95e4-8978dbb3d590.root",
            path + "/DstarToD0Pi_D0To2Mu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/4264fdbb-4822-43bb-966b-2a3da02bc49b.root",
            path + "/DstarToD0Pi_D0To2Mu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/5c7f65e6-9a42-4e87-9463-32b437dfd863.root",
            path + "/DstarToD0Pi_D0To2Mu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/7ff86924-1e30-4633-ba52-0915badd3f66.root",
            path + "/DstarToD0Pi_D0To2Mu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/95982060-e6de-4d62-b0e2-15dd9adaac06.root",
            path + "/DstarToD0Pi_D0To2Mu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/a121d5c0-d816-42e6-a415-31706a143da1.root",
            path + "/DstarToD0Pi_D0To2Mu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/ae31daa2-7105-4bdc-b456-c1c10d9a81ba.root",
            path + "/DstarToD0Pi_D0To2Mu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/af37d534-abfe-4ad8-935a-22c6971b6ef8.root",
            path + "/DstarToD0Pi_D0To2Mu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/b9c0f04e-75d8-489a-9670-b1e720ce532e.root",
            path + "/DstarToD0Pi_D0To2Mu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/bccd752e-2f56-494b-83f0-86fd0d8059c4.root",
            path + "/DstarToD0Pi_D0To2Mu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/c7ab600c-b132-42d1-a07c-4dab362332bc.root",
            path + "/DstarToD0Pi_D0To2Mu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/d0400b0b-f754-42ee-ba1e-cd9838c2a4aa.root",
            path + "/DstarToD0Pi_D0To2Mu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/d3b0ec5b-97c6-4b66-8370-c61de101327a.root",
            path + "/DstarToD0Pi_D0To2Mu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/dbf39517-9639-4f68-9fb9-eecf8c21e431.root",
            path + "/DstarToD0Pi_D0To2Mu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/de11d9f6-f03f-4b76-806c-65120b70cf36.root",
            path + "/DstarToD0Pi_D0To2Mu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/e86f4986-bd7d-4739-acd9-c4da17b471ec.root",
            path + "/DstarToD0Pi_D0To2Mu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/ea54cd31-9792-4d65-8c21-1236921055ff.root",
            path + "/DstarToD0Pi_D0To2Mu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/eb9dfbdc-addd-4c63-903d-8edc27c33428.root",
            path + "/DstarToD0Pi_D0To2Mu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/f5217816-f36b-441a-8d28-6112aeab01a7.root",
        ],
        'chain':None
    },
    {
        'final_state':'pipi',
        'name':'\dzpipi',
        'files':[
            "/eos/cms/store/group/phys_muon/dmytro/tmp/DstarToD0Pi_D0To2Pi_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen.root",
        ],
        'chain':None
    },
]

cuts = [
    {
        'cut':{
            'mm':'dstar_gen_pdgId!=0 && dstar_mm_index>=0 && mm_mu1_pt[dstar_mm_index]>4 && mm_mu2_pt[dstar_mm_index]>4',
            'pipi':'dstar_hh_index>=0&&abs(dstar_gen_pdgId)==413&&abs(hh_gen_had1_pdgId)==211&&abs(hh_gen_had2_pdgId)==211 && hh_had1_pt[dstar_hh_index]>4 && hh_had2_pt[dstar_hh_index]>4',
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
            'fakes':'dstar_dm_pv>0.140 && dstar_dm_pv<0.155',
        },
        'name':'$\dm$ in [0.140, 0.155]',
    },
    {
        'cut':{
            'mm':'mm_kin_mass[dstar_mm_index]>1.81 && mm_kin_mass[dstar_mm_index]<1.94',
            'pipi':'hh_kin_mass[dstar_hh_index]>1.81 && hh_kin_mass[dstar_hh_index]<1.94',
            'fakes':'mm_kin_mass[dstar_mm_index]>1.81 && mm_kin_mass[dstar_mm_index]<1.94',
        },
        'name':'$\PDz$ mass in [1.81, 1.94]',
    },
    {
        'cut':{
            'mm':'mm_kin_vtx_prob[dstar_mm_index]>0.01',
            'pipi':'hh_kin_vtx_prob[dstar_hh_index]>0.01',
            'fakes':'mm_kin_vtx_prob[dstar_mm_index]>0.01',
        },
        'name':'$\PDz$ vertex probability $> 0.01$',
    },
    {
        'cut':{
            'mm':'dstar_pv_with_pion_prob>0.1',
            'pipi':'dstar_pv_with_pion_prob>0.1',
            'fakes':'dstar_pv_with_pion_prob>0.1',
        },
        'name':'$\PDstpm$ vertex probability $> 0.1$',
    },
    {
        'cut':{
            'mm':'Muon_softMva[mm_mu1_index[dstar_mm_index]] > 0.45 && Muon_softMva[mm_mu2_index[dstar_mm_index]] > 0.45',
            'fakes':'Muon_softMva[mm_mu1_index[dstar_mm_index]] > 0.45 && Muon_softMva[mm_mu2_index[dstar_mm_index]] > 0.45',
        },
        'name':'Muon identification',
    },
    {
        'cut':{
            'mm':'mm_kin_sl3d[dstar_mm_index]>3',
            'pipi':'hh_kin_sl3d[dstar_hh_index]>3',
            'fakes':'mm_kin_sl3d[dstar_mm_index]>3',
        },
        'name':r'$\fls > 3$',
    },
    {
        'cut':{
            'mm':'mm_kin_alpha[dstar_mm_index]<0.1',
            'pipi':'hh_kin_alpha[dstar_hh_index]<0.1',
            'fakes':'mm_kin_alpha[dstar_mm_index]<0.1',
        },
        'name':'$\pa < 0.1$',
    },
    
]

report = EfficiencyReport(samples, cuts)
report.make_report("gen", r"%7.4f")

