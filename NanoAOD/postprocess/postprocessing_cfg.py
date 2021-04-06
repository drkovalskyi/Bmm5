# Config file for the postprocessor

workdir = "/afs/cern.ch/work/d/dmytro/projects/RunII-NanoAODv6/src/Bmm5/NanoAOD/postprocess/"
version = 511
input_location = "/eos/cms/store/group/phys_bphys/bmm/bmm5/NanoAOD"
output_location = "/eos/cms/store/group/phys_bphys/bmm/bmm5/PostProcessing"
xrootd_prefix = "root://eoscms.cern.ch:/"
web_report_path = "/afs/cern.ch/user/d/dmytro/www/public_html/bmm5/postprocessing/"

odebug = False
tmp_prefix = "tmpPPNA"
require_log_for_success = True

# Test RegEx at https://www.regextester.com/
tasks = [
    {
        'input_pattern':'EGamma|DoubleEG|SingleElectron|SinglePhoton|MINIAODSIM',
        # 'input_pattern':'EGamma|DoubleEG|SingleElectron|SinglePhoton',
        'processor':'Skimmer',
        'cut':'nks>0',
        'name':'ks',
        'type':'NanoAOD-skims',
        'files_per_job':100
    },
    {
        'input_pattern':'EGamma|DoubleEG|SingleElectron|SinglePhoton|MINIAODSIM',
        # 'input_pattern':'EGamma|DoubleEG|SingleElectron|SinglePhoton',
        'processor':'Skimmer',
        'cut':'nlambda>0',
        'name':'lambda',
        'type':'NanoAOD-skims',
        'files_per_job':100
    },
    {
        'input_pattern':'EGamma|DoubleEG|SingleElectron|SinglePhoton|MINIAODSIM',
        # 'input_pattern':'EGamma|DoubleEG|SingleElectron|SinglePhoton',
        'processor':'Skimmer',
        'cut':'((phi_trk1_pt > 4 && phi_trk2_pt > 3.0) || ( phi_trk2_pt > 4 && phi_trk1_pt > 3.0)) &&  phi_kin_vtx_prob>0.3 && phi_kin_sipPV<1  && phi_doca<0.004 && phi_kin_lxy<4',
        'name':'phi',
        'type':'NanoAOD-skims',
        'files_per_job':100
    },
    {
        'input_pattern':'EGamma|DoubleEG|SingleElectron|SinglePhoton|MINIAODSIM',
        # 'input_pattern':'EGamma|DoubleEG|SingleElectron|SinglePhoton',
        'processor':'Skimmer',
        'cut':'phi_ds_pion_pt>0',
        'name':'ds',
        'type':'NanoAOD-skims',
        'files_per_job':100
    },
    {
        'input_pattern':'EGamma|DoubleEG|SingleElectron|SinglePhoton',
        'processor':'Skimmer',
        'cut':'nmm>0',
        'name':'mm',
        'type':'NanoAOD-skims',
        'files_per_job':100
    },
    {
        'input_pattern':'Charmonium',
        'processor':'FlatNtupleForBmmMva',
        'signal_only': False,
        'tree_name': "mva",
        'blind': True,
        'name':'bmm_mva',
        'type':'FlatNtuples',
        'files_per_job':20
    },
    {
        'input_pattern':'BsToMuMu_BMuonFilter',
        'processor':'FlatNtupleForBmmMva',
        'signal_only': True,
        'tree_name': "mva",
        'blind': False,
        'name':'bmm_mva',
        'type':'FlatNtuples',
        'files_per_job':20
    },
]

resources = [
    "LocalResourceHandler(16)", #vocms0118.cern.ch
    "SSHResourceHandler('vocms0109.cern.ch',32)",
]
