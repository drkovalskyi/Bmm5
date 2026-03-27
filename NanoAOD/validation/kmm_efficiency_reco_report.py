import ROOT
import glob

from base_efficiency_reco_report import EfficiencyReport

limit = 1000
path = "/eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/532"
path2 = "/eos/cms/store/group/phys_bphys/bmm/bmm6/PostProcessing/Skims/531/"
path3 = "/eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/531"
samples = [
    {
        'final_state':'mm',
        'name':'\ksmm',
        'scale':1e4,
        # 'files': glob.glob(f'{path}/K0sToMuMu_K0sFilter_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v14-v4+MINIAODSIM/*.root')[:limit],
        'files': glob.glob(f'{path}/K0sToMuMu_K0sFilter_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v1+MINIAODSIM/*.root')[:limit],
        'chain':None
    },
    {
        'final_state':'pipi',
        'name':'\kspipi',
        'scale':1e4,
        # 'scale':302, # special value to match Ksmm reco efficiency
        # 'files': glob.glob(f'{path2}/kspipi_gen/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8+Run3Summer22MiniAODv4-validDigi_130X_mcRun3_2022_realistic_v5-v4+MINIAODSIM/*.root')[:limit],
        'files': glob.glob(f'{path}/K0sToPiPi_K0sFilter_TuneCP5_13p6TeV_pythia8-evtgen/*.root')[:limit],
        'chain':None
    },
    {
        'final_state':'fakes',
        'name':r'\kspipimmshort',
        'scale':1e10/50/50,
        # 'files': glob.glob(f'{path1}/K0sTo2PiTo2Mu_K0sFilter_MuFilter_PiLifetime0p02_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v14-v2+MINIAODSIM/*.root')[:limit],
        'files': glob.glob(f'{path3}/K0sTo2PiTo2Mu_K0sFilter_MuFilter_PiLifetime0p02_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v2+MINIAODSIM/*.root')[:limit],
        'chain':None
    },
]

cuts = [
    {
        'cut':{
            'mm':'abs(mm_gen_ks_pdgId)==310 && mm_mu1_pt>4 && mm_mu2_pt>3',
            'pipi':'abs(hh_gen_ks_pdgId)==310 && hh_had1_pt>4 && hh_had2_pt>3',
            'fakes':'abs(mm_gen_cpdgId)==310 && mm_mu1_pt>4 && mm_mu2_pt>3',
        },
        'name':'MC matched reco with quality requirements',
    },
    {
        'cut':{
            'mm':'mm_kin_vtx_prob>0.01',
            'pipi':'hh_kin_vtx_prob>0.01',
            'fakes':'mm_kin_vtx_prob>0.01'
        },
        'name':'Vertex probability $>0.01$',
    },
    {
        'cut':{
            'mm':'mm_kin_slxy > 3',
            'pipi':'hh_kin_slxy > 3',
            'fakes':'mm_kin_slxy > 3'
        },
        'name':'Tranverse vertex displacement significance $>3$',
    },
    {
        'cut':{
            'mm':'mm_kin_lxy > 1',
            'pipi':'hh_kin_lxy > 1',
            'fakes':'mm_kin_lxy > 1'
        },
        'name':'Tranverse vertex displacement $>1$ cm',
    },
    {
        'cut':{
            'mm':'mm_kin_alpha<0.1',
            'pipi':'hh_kin_alpha<0.1',
            'fakes':'mm_kin_alpha<0.1',
        },
        'name':'Pointing agnle $< 0.1$',
    },
    {
        'cut':{
            'mm':'Muon_mediumId[mm_mu1_index] && Muon_mediumId[mm_mu2_index]',
            'fakes':'Muon_mediumId[mm_mu1_index] && Muon_mediumId[mm_mu2_index]',
        },
        'name':'Medium Muon ID',
    },
    {
        'cut':{
            'mm':'HLT_DoubleMu4_3_LowMass',
            'fakes':'HLT_DoubleMu4_3_LowMass',
        },
        'name':r'\verb|HLT_DoubleMu4_3_LowMass|',
    },
    
]

report = EfficiencyReport(samples, cuts)
report.make_report("gen", r"%7.4f")

