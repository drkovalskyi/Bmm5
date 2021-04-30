from ModelHandler import *
from train_muon_mva import MuonMVA

fig = None
ax = None

def setup():
    global fig
    global ax
    fig, ax = plt.subplots()
    # help(fig)
    # help(ax)
    plt.grid(which='both', axis='both')

def save_result(name):
    ax.set_xlabel("bkg eff")
    ax.set_ylabel("sig eff")
    ax.set_xlim([0.0, 1.0])
    ax.set_ylim([0.7, 1.0])
    ax.set_title("ROC curves")
    ax.legend(loc='lower right')
    fig.set_tight_layout(True)
    # ax.set_yscale("log")
    # ax.set_xscale("log")
    fig.savefig("%s.pdf"%name)

def load_model(files, model_name, n_split=3, event_index=None):
    if event_index != None:
        model_name = "%s-Event%u" % (model_name, event_index)

    features = json.load(open("%s.features" % (model_name)))
    model = MuonMVA("",features, "muons", "evt")
    model.load_datasets(files)
    model.apply_selection(model.data["pt"]>5)
    model.apply_selection(model.data["pt"]<6)
    model.apply_selection(abs(model.data["eta"])<1.4)

    model.bst = xgb.Booster({'nthread': 4})  # init model
    model.bst.load_model("%s.model" % (model_name))
    
    if event_index!=None:
        model.prepare_train_and_test_datasets(n_split, event_index)
    else:
        model.y_train = model.get_true_classification(model.data)
        model.y_test  = model.get_true_classification(model.data)
        model.x_train = model.get_feature_data(model.data)
        model.x_test  = model.get_feature_data(model.data)
    return model
    
def add_roc_curve(model, label):
    dtest = xgb.DMatrix(model.x_test, label=model.y_test, feature_names=model.features)
    y_pred_test = model.bst.predict(dtest)

    fpr_test, tpr_test, _ = roc_curve(model.y_test, y_pred_test)
    ax.plot(fpr_test, tpr_test, label=label)

def plot_efficiency_pt(model, label):
    dtest = xgb.DMatrix(model.x_test, label=model.y_test, feature_names=model.features)
    y_pred_test = model.bst.predict(dtest)

    fpr_test, tpr_test, _ = roc_curve(model.y_test, y_pred_test)
    ax.plot(fpr_test, tpr_test, label=label)

    
def add_old_soft_mva_roc_curve_without_preselection(model):
    data = model.test_data
    if not data:
        data = model.data
    y_bdt = (data["softMva"] )
    fpr_ref, tpr_ref, _ = roc_curve(model.y_test , y_bdt)
    ax.plot(fpr_ref, tpr_ref, label="Bmm4 2016 soft Muon MVA")

def get_sig_eff(model, bkg_eff):
    dtest = xgb.DMatrix(model.x_test, label=model.y_test, feature_names=model.features)
    y_pred_test = model.bst.predict(dtest)

    sig = []
    bkg = []
    for i, y in enumerate(model.y_test):
        if y:
            sig.append(y_pred_test[i])
        else:
            bkg.append(y_pred_test[i])
    
    quantile = np.percentile(bkg, (1-bkg_eff)*100.)
    print "quantile:", quantile
    n_sig = 0
    for entry in sig:
        if entry >= quantile:
            n_sig += 1
    return float(n_sig)/len(sig)

def add_standard_selectors(model):
    vars = {'softMvaId':"*", "mediumId":"o"}
    # print data
    data = model.test_data
    if not data:
        data = model.data
    for var, marker in vars.items():
        print var
        y_data = model.select_events(data, data["sim_type"] == 1)
        x_data = model.select_events(data, data["sim_type"] != 1)
        y_passed = model.select_events(y_data, y_data[var]==1)
        x_passed = model.select_events(x_data, x_data[var]==1)
        bkg_eff = float(len(y_passed["sim_type"]))/len(y_data["sim_type"])
        print "bkg_eff:", bkg_eff
        sig_eff = float(len(x_passed["sim_type"]))/len(x_data["sim_type"])
        print "sig_eff:", sig_eff
        model_sig_eff = get_sig_eff(model, bkg_eff)
        print "model sig_eff:", model_sig_eff

        ax.plot([bkg_eff], [sig_eff], marker=marker, label=var)

files = [
    "/eos/cms/store/group/phys_bphys/bmm/bmm5/PostProcessing/FlatNtuples/512/muon_mva/QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_13TeV_pythia8+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v3+MINIAODSIM/15475f7d44f873b2f0aaca4372b12284.root"
    # "/eos/cms/store/group/phys_bphys/bmm/bmm5/PostProcessing/FlatNtuples/512/muon_mva/QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_13TeV_pythia8+RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1+MINIAODSIM/0ebc8af01e485c5e733523cc21138630.root"
    # "/eos/cms/store/group/phys_bphys/bmm/bmm5/PostProcessing/FlatNtuples/512/muon_mva/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8+RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2+MINIAODSIM/4969bbd7ff55a80c1d740a89c64a6b60.root"
]

setup()

# model = load_model(files, "results/muon_mva/Run2018-20210421-0013", 3, 0)
# # model.apply_selection(model.data["pt"]>5)
# add_roc_curve(model, "Run2018-20210421-0013")
# add_old_soft_mva_roc_curve_without_preselection(model)
# add_standard_selectors(model)

# model = load_model(files, "results/muon_mva/Run2017-20210422-1123", 3, 0)
# add_roc_curve(model, "Run2017-20210422-1123")
# add_old_soft_mva_roc_curve_without_preselection(model)
# add_standard_selectors(model)

model = load_model(files, "results/muon_mva/Run2018-20210430-0757", 3, 0)
add_roc_curve(model, "Run2018-20210430-0757")
# model = load_model(files, "results/muon_mva/Run2018-20210424-1707", 3, 0)
# add_roc_curve(model, "Run2018-20210424-1707")
# model = load_model(files, "results/muon_mva/Run2017-20210422-1123", 3, 0)
# add_roc_curve(model, "Run2017-20210422-1123")
# model = load_model(files, "results/muon_mva/Run2016-20210422-1123", 3, 0)
# add_roc_curve(model, "Run2016-20210422-1123")

add_old_soft_mva_roc_curve_without_preselection(model)
add_standard_selectors(model)

# model = load_model(files, "results/muon_mva/Run2016-20210422-1123", 3, 0)
# add_roc_curve(model, "Run2016-20210422-1123")
# add_old_soft_mva_roc_curve_without_preselection(model)
# add_standard_selectors(model)




# model = load_model(files, "results/muon_mva/Run2018-20210422-0153", 3, 0)
# add_roc_curve(model, "Run2018-20210422-0153")

save_result("performance_muon_mva_bkg_vs_sig")
