from ModelHandler import *

feature_names = [
    #
    # Features in order of decreasing contribution
    #
    ### the best
    "mm_kin_alpha",        # used in old BDT
    "mm_kin_alphaBS",
    "mm_kin_spvip",        # used in old BDT
    "mm_kin_pvip",
    "mm_iso",              # used in old BDT
    ### good
    "mm_m1iso",            # used in old BDT
    "mm_m2iso",            # used in old BDT
    "mm_kin_sl3d",         # used in old BDT
    "mm_kin_vtx_chi2dof",  # used in old BDT
    ### weak
    "mm_nBMTrks",
    # "mm_mu2_pt",
    # "mm_closetrks1",
    # "mm_nDisTrks",
    # "mm_kin_pt"
    ### very weak (in order of decreasing contribution)
    # "mm_mu2_eta",
    # "mm_docatrk",     # used in old BDT
    # "mm_closetrk",    # used in old BDT
    # "mm_kin_eta",     # used in old BDT
    # "mm_kin_spvlip",
    # "mm_kin_pvlip",
    # "mm_otherVtxMaxProb",
    ### new
    "mm_otherVtxMaxProb1",
    "mm_otherVtxMaxProb2",

]

class BmmMVA(ModelHandler):
    """ Bmm MVA specific parameters and methods """
    def __init__(self, output_path, features, tree_name, event_branch_name):
        ModelHandler.__init__(self, output_path, features, tree_name, event_branch_name)
        self.roc_bkg_min = 1e-6
        self.roc_bkg_max = 1.0
        self.roc_sig_min = 0.1
        self.roc_sig_max = 1.0
        
    def get_true_classification(self, data):
        """ Define signal selection"""
        return (data["mm_gen_pdgId"] != 0)

    def get_parameters(self):
        parameters = super(BmmMVA, self).get_parameters()

        ## Fast
        # parameters['min_child_weight'] = 0.01
        # parameters['max_depth'] = 3
        # parameters['eta'] = 1

        ## Best
        parameters['min_child_weight'] = 0.00001
        parameters['max_depth'] = 10
        parameters['eta'] = 0.01
        # parameters['nthread'] = 30 

        parameters['tree_method'] = 'hist' 
        return parameters

def train_model(files, name):
    
    model = BmmMVA("results/bmm_mva", feature_names, "mva", "evt_event")
    model.load_datasets(files)
    model.apply_selection(model.data["HLT_DoubleMu4_3_Bs"]==1)

    number_of_splits = 3
    train_model_for_all_splits = True
    for event_index in range(number_of_splits):
        model.prepare_train_and_test_datasets(number_of_splits, event_index)
        print "Number of signal/background events in training sample: %u/%u (%0.3f)" \
            % (sum(model.y_train == True),
               sum(model.y_train == False),
               sum(model.y_train == True) / float(sum(model.y_train == False)))
        model_name = "%s-%s-Event%u" % (name, date, event_index)
        model.train(model_name)
        if not train_model_for_all_splits:
            break

if __name__ == "__main__":
    date = datetime.now().strftime("%Y%m%d-%H%M")

    data_path = "/eos/cms/store/group/phys_bphys/bmm/bmm5/PostProcessing/FlatNtuples/516/bmm_mva/"
    
    files_Run2018 = [
        data_path + "BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v1+MINIAODSIM/",
        data_path + "Charmonium+Run2018A-12Nov2019_UL2018_rsb-v1+MINIAOD/",
        data_path + "Charmonium+Run2018B-12Nov2019_UL2018-v1+MINIAOD/",
        data_path + "Charmonium+Run2018C-12Nov2019_UL2018_rsb_v3-v1+MINIAOD/",
        data_path + "Charmonium+Run2018D-12Nov2019_UL2018-v1+MINIAOD/",
    ]
    files_Run2017 = [
        data_path + "BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL17MiniAOD-106X_mc2017_realistic_v6-v1+MINIAODSIM/",
        data_path + "Charmonium+Run2017B-09Aug2019_UL2017-v1+MINIAOD/",
        data_path + "Charmonium+Run2017C-09Aug2019_UL2017-v1+MINIAOD/",
        data_path + "Charmonium+Run2017D-09Aug2019_UL2017-v1+MINIAOD/",
        data_path + "Charmonium+Run2017E-09Aug2019_UL2017-v1+MINIAOD/",
        data_path + "Charmonium+Run2017F-09Aug2019_UL2017-v1+MINIAOD/",
    ]

    # train_model(files_Run2018, "Run2018")

    train_model(files_Run2017 + files_Run2018, "Run2017-2018")

    
    # files_Run2017 = [
    #     data_path + "QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_13TeV_pythia8+RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1+MINIAODSIM/",
    #     # data_path + "QCD_Pt-50to80_MuEnrichedPt5_TuneCP5_13TeV_pythia8+RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1+MINIAODSIM/",
    #     # data_path + "QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_13TeV_pythia8+RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1+MINIAODSIM/"
    # ]
    # train_model(files_Run2017, "Run2017")

    # files_Run2016 = [
    #     data_path + "QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8+RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2+MINIAODSIM/",
    #     # data_path + "QCD_Pt-50to80_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8+RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2+MINIAODSIM/",
    #     # data_path + "QCD_Pt-80to120_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8+RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2+MINIAODSIM/",
    #     # data_path + "QCD_Pt-80to120_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8+RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2+MINIAODSIM/"
    # ]
    # train_model(files_Run2016, "Run2016")
    

