# import os
import sys, pickle, json
from datetime import datetime

import numpy as np
import uproot
from sklearn.model_selection import train_test_split

# from tqdm import tqdm
import xgboost as xgb
from sklearn.metrics import roc_curve, roc_auc_score

from pprint import pprint

import matplotlib as mpl
# https://matplotlib.org/faq/usage_faq.html
mpl.use('Agg')
import matplotlib.pyplot as plt

output_path = "results/"

def load_data(fname):
    """Load ROOT data and turn MVA tree into an array.

    The array is represented in a form of dictionary with the branch
    name representing a key and the value is numpy.ndarray. 
    """
    print "Loading data from", fname
    f = uproot.open(fname)
    tree = f["mva"]
    return tree.arrays(tree.keys())

def merge_data(datasets):
    """Merge multiple datasets into one dataset."""
    dataset = dict()
    if len(datasets)<1: return dataset
    branches = datasets[0].keys()
    for branch in branches:
        dataset[branch] = np.concatenate([ds[branch] for ds in datasets],axis=None)
    return dataset

def load_datasets(file_list=[]):
    """Load and merge ROOT trees with MVA data into a single dataset."""
    datasets = []
    for f in file_list:
        datasets.append(load_data(f))
    return merge_data(datasets)


def get_true_classification(data):
    """Make a column representing the true classification of events.
    
    The collumn is used for training and performance testing.
    """
    y_data = (data["mm_gen_pdgId"]!=0 )
    return y_data

def get_feature_data(data,feature_names):
    """Get columns containing feature data."""
    x_data = np.column_stack([data[name] for name in feature_names])
    return x_data

def get_train_and_test_datasets_split_randomly(data):
    """Split dataset into random train and test subsets."""
    y_data = get_true_classification(data)
    x_data = get_feature_data(data, feature_names)
    x_train, x_test, y_train, y_test = train_test_split(x_data, y_data, test_size=0.25, random_state=42)
    return x_train,x_test,y_train,y_test

def select_events(data,condition):
    """Select events"""
    selected_data = dict()
    indices = np.where(condition)[0]
    for branch in data:
        selected_data[branch] = data[branch].take(indices)
    return selected_data

def get_train_and_test_datasets_split_be_event_number(data,n=4,index=0):
    """Split dataset into train and test subsets using event number.
    
    When training is performed on data used in the analysis for other
    purposes, one needs to make sure that no bias is introduced by
    using the same events for training and classification. This can be
    achieved by performing multiple training splitting events in
    non-overlapping groups based on the modulus of the event number
    with the divisor n. By default n=4. This means that one needs to
    perform 4 trainings with different index number. In each training
    (1-1/4)=75% of available data with modulus not equal to index is
    used in the training data set.
    """
    train_data = select_events(data,data['evt_event']%n!=index)
    test_data  = select_events(data,data['evt_event']%n==index)
    y_train = get_true_classification(train_data)
    y_test  = get_true_classification(test_data)
    x_train = get_feature_data(train_data, feature_names)
    x_test  = get_feature_data(test_data,  feature_names)
    
    return x_train,x_test,y_train,y_test

def get_parametrs():
    """Get parameters for XGBoost."""
    param = {}
    # https://xgboost.readthedocs.io/en/latest/parameter.html
    # objective [default=reg:squarederror]
    #
    # - binary:logistic: logistic regression for binary classification,
    # output probability
    param['objective'] = 'binary:logistic'
    # eta [default=0.3, alias: learning_rate] 
    #
    # - Step size shrinkage used in update to prevents overfitting. After
    # each boosting step, we can directly get the weights of new features,
    # and eta shrinks the feature weights to make the boosting process
    # more conservative.
    # - range: [0,1]
    param['eta'] = 0.01
    # max_depth [default=6]
    #
    # - Maximum depth of a tree. Increasing this value will make the model
    # more complex and more likely to overfit. 0 is only accepted in
    # lossguided growing policy when tree_method is set as hist and it
    # indicates no limit on depth. Beware that XGBoost aggressively
    # consumes memory when training a deep tree.  
    # - range: [0,Inf] (0 is only accepted in lossguided growing policy
    # when tree_method is set as hist)
    param['max_depth'] = 10
    param['silent'] = 1
    # number of threads to use for training. It's wise to set it to
    # the number of real CPUs. More is not always better.
    param['nthread'] = 15
    # eval_metric [default according to objective]
    # 
    # - Evaluation metrics for validation data, a default metric will be
    # assigned according to objective (rmse for regression, and error for
    # classification, mean average precision for ranking)
    # - auc: Area under the curve
    param['eval_metric'] = "auc"
    # param['eval_metric'] = "aucpr"
    # subsample [default=1]
    # 
    # - Subsample ratio of the training instances. Setting it to 0.5 means
    # that XGBoost would randomly sample half of the training data prior
    # to growing trees. and this will prevent overfitting. Subsampling
    # will occur once in every boosting iteration.
    # - range: (0,1]
    param['subsample'] = 0.6
    # alpha [default=0, alias: reg_alpha]
    #
    # - L1 regularization term on weights. Increasing this value will make
    # model more conservative. Normalised to number of training examples.
    param['alpha'] = 8.0
    # gamma [default=0, alias: min_split_loss]
    # 
    # - Minimum loss reduction required to make a further partition on a
    # leaf node of the tree. The larger gamma is, the more conservative
    # the algorithm will be.
    # - range: [0,inf]
    param['gamma'] = 2.0
    # lambda [default=0, alias: reg_lambda]
    # 
    # - L2 regularization term on weights. Increasing this value will make
    # model more conservative. Normalised to number of training examples.
    param['lambda'] = 1.0
    # min_child_weight [default=1]
    #
    # - Minimum sum of instance weight (hessian) needed in a child. If the
    # tree partition step results in a leaf node with the sum of instance
    # weight less than min_child_weight, then the building process will
    # give up further partitioning. In linear regression task, this simply
    # corresponds to minimum number of instances needed to be in each
    # node. The larger min_child_weight is, the more conservative the
    # algorithm will be.
    # - range: [0,Inf]
    # param['min_child_weight'] = 0.01
    param['min_child_weight'] = 0.00001
    # colsample_bytree, colsample_bylevel, colsample_bynode [default=1]
    # 
    # - This is a family of parameters for subsampling of columns.
    # - All colsample_by* parameters have a range of (0, 1], the default
    # value of 1, and specify the fraction of columns to be subsampled.
    # - colsample_bytree is the subsample ratio of columns when constructing
    # each tree. Subsampling occurs once for every tree constructed.
    param['colsample_bytree'] = 1.0
    
    return param

def train(x_train,x_test,y_train,y_test,model_name="test"):
    """Train a model."""
    dtrain = xgb.DMatrix( x_train, label=y_train, feature_names=feature_names)
    dtest = xgb.DMatrix( x_test, label=y_test, feature_names=feature_names)
    evallist  = [(dtrain,'train'), (dtest,'eval')]

    num_round = 5000

    param = get_parametrs()

    sumw_pos = np.abs(dtrain.get_label()==1).sum()
    sumw_neg = np.abs(dtrain.get_label()==0).sum()
    # scale_pos_weight [default=1]
    #
    # - Control the balance of positive and negative weights, useful for
    # unbalanced classes. A typical value to consider: sum(negative
    # instances) / sum(positive instances). See Parameters Tuning for more
    # discussion.
    # https://xgboost.readthedocs.io/en/latest/tutorials/param_tuning.html

    # target_s_over_b = 0.005
    target_s_over_b = .1

    param["scale_pos_weight"] = sumw_neg/float(sumw_pos)*target_s_over_b
    pprint(param)

    bst = xgb.train( param.items(), dtrain, num_round, evallist, early_stopping_rounds=10 )
    # bst.feature_names = feature_names
    max_length = max([len(s) for s in feature_names])
    feature_format = "%"+"%u"%max_length+"s"
    scores = bst.get_score(importance_type='gain')
    print "Importance scores:"
    for name in feature_names:
        if name in scores:
            print (feature_format+" %0.1f") % (name,scores[name])
        else:
            print (feature_format+" unused") % name
    pickle.dump({'model':bst,'feature_names':feature_names},open("%s/%s.pkl"%(output_path,model_name),"wb"))
    bst.dump_model("%s/%s.txt"%(output_path,model_name))
    bst.save_model("%s/%s.model"%(output_path,model_name))
    json.dump(feature_names,open("%s/%s.features"%(output_path,model_name),"w"))
    json.dump(param,open("%s/%s.params"%(output_path,model_name),"w"))

    # write_json("model.json", bst, feature_names)
    # json_to_cfunc("model.json", fname_out="func.h")
    return bst,dtrain,dtest

def cross_validate(x_train,y_train,model_name="test"):
    """Cross-validate model and its parameters"""
    dtrain = xgb.DMatrix( x_train, label=y_train, feature_names=feature_names)

    num_round = 200

    param = get_parametrs()

    sumw_pos = np.abs(dtrain.get_label()==1).sum()
    sumw_neg = np.abs(dtrain.get_label()==0).sum()
    # scale_pos_weight [default=1]
    #
    # - Control the balance of positive and negative weights, useful for
    # unbalanced classes. A typical value to consider: sum(negative
    # instances) / sum(positive instances). See Parameters Tuning for more
    # discussion.
    # https://xgboost.readthedocs.io/en/latest/tutorials/param_tuning.html

    target_s_over_b = 0.003

    param["scale_pos_weight"] = sumw_neg/float(sumw_pos)*target_s_over_b
    pprint(param)

    cv_results = xgb.cv( 
        param.items(), 
        dtrain, 
        num_round, 
        seed=42, # fix seed to be able to compare scores with different parameters
        nfold=5, # the number of folds to use for cross-validation
        metrics={'auc'},
        verbose_eval=True,
        early_stopping_rounds=10 )
    
    pprint(cv_results)
    # # bst.feature_names = feature_names
    # max_length = max([len(s) for s in feature_names])
    # feature_format = "%"+"%u"%max_length+"s"
    # scores = bst.get_score(importance_type='gain')
    # print "Importance scores:"
    # for name in feature_names:
    #     if name in scores:
    #         print (feature_format+" %0.1f") % (name,scores[name])
    #     else:
    #         print (feature_format+" unused") % name
    # pickle.dump({'model':bst,'feature_names':feature_names},open("%s.pkl"%model_name,"wb"))
    # bst.dump_model("%s.txt"%model_name)
    # bst.save_model("%s.model"%model_name)
    # json.dump(feature_names,open("%s.features"%model_name,"w"))

    # # write_json("model.json", bst, feature_names)
    # # json_to_cfunc("model.json", fname_out="func.h")
    # return bst,dtrain,dtest


def make_validation_plots(bst,dtrain,dtest,model_name="test"):
    y_pred_train = bst.predict(dtrain)
    y_pred_test  = bst.predict( dtest)

    fig, ax = plt.subplots()
    preds_bkg_test  = y_pred_test[y_test==0]
    preds_sig_test  = y_pred_test[y_test==1]
    preds_bkg_train = y_pred_train[y_train==0]
    preds_sig_train = y_pred_train[y_train==1]
    
    density = True
    bins = np.linspace(0.0,1,50)
    edges = 0.5*(bins[:-1]+bins[1:])
    
    from scipy import stats
    pks_bkg = stats.ks_2samp(preds_bkg_train,preds_bkg_test)[1]
    pks_sig = stats.ks_2samp(preds_sig_train,preds_sig_test)[1]

    counts_train_bkg,_,_ = ax.hist(preds_bkg_train, bins=bins,histtype="stepfilled",alpha=0.45, normed=density, label="bkg, train",color=["b"])
    counts_train_sig,_,_ = ax.hist(preds_sig_train, bins=bins,histtype="stepfilled",alpha=0.45, normed=density, label="sig, train",color=["r"])
    counts_test_bkg,_,_ = ax.hist(preds_bkg_test, bins=bins,histtype="step",alpha=1.0, normed=density, label="bkg, test (KS prob = {:.2f})".format(pks_bkg),color=["b"], lw=1.5, linestyle="solid")
    counts_test_sig,_,_= ax.hist(preds_sig_test, bins=bins,histtype="step",alpha=1.0, normed=density, label="sig, test (KS prob = {:.2f})".format(pks_sig),color=["r"], lw=1.5, linestyle="solid")

    ax.set_yscale("log")
    ax.legend()

    ax.set_ylim([0.01,ax.get_ylim()[1]])
    fig.set_tight_layout(True)
    fig.savefig("%s/%s-validation-disc.pdf"%(output_path,model_name))

    fpr_test,tpr_test,_ = roc_curve(y_test, y_pred_test)
    fpr_train,tpr_train,_ = roc_curve(y_train, y_pred_train)
    auc_test = roc_auc_score(y_test, y_pred_test)
    auc_train = roc_auc_score(y_train, y_pred_train)
    fig, ax = plt.subplots()
    plt.grid(which='both', axis='both')
    ax.plot(fpr_test,tpr_test, label="test AUC = {:.6f}".format(auc_test))
    ax.plot(fpr_train,tpr_train, label="train AUC = {:.6f}".format(auc_train))
    ax.set_xlabel("bkg eff")
    ax.set_ylabel("sig eff")
    ax.set_xlim([0.00001,0.1])
    ax.set_ylim([0.2,1.0])
    ax.set_title("ROC curves")
    ax.legend()
    fig.set_tight_layout(True)
    ax.set_yscale("log")
    ax.set_xscale("log")
    fig.savefig("%s/%s-validation-roc.pdf"%(output_path,model_name))

def prepare_dataset():
    all_data = load_datasets([
        # "python/BdToMuMu_BMuonFilter_RunIIAutumn18MiniAOD.root",
        # "python/BsToMuMu.root",
        "python/BsToMuMu_BMuonFilter_RunIIAutumn18MiniAOD.root",
        "python/Charmonium+Run2018A.root",
        "python/Charmonium+Run2018B.root",
        "python/Charmonium+Run2018C.root",
        "python/Charmonium+Run2018D.root",
        # "python/BsToMuMu_BMuonFilter_RunIIFall17MiniAODv2.root",
        "python/Charmonium+Run2017B.root",
        "python/Charmonium+Run2017C.root",
        "python/Charmonium+Run2017D.root",
        "python/Charmonium+Run2017E.root",
        "python/Charmonium+Run2017F.root",
        # "python/QCD_Pt-30to50_MuEnrichedPt5_RunIIAutumn18MiniAOD.root",
        # "python/QCD_Pt-50to80_MuEnrichedPt5_RunIIAutumn18MiniAOD.root",
        # "python/QCD_Pt-80to120_MuEnrichedPt5_RunIIAutumn18MiniAOD.root",
    ])
    print "Total number of events:", len(all_data['evt_event'])

    data = select_events(all_data,all_data['HLT_DoubleMu4_3_Bs']>0)

    print "Number of events selected for training:", len(data['evt_event'])
        
    return data
    
feature_names = [
    #
    # Features in order of decreasing contribution
    #
    ### the best
    "mm_kin_alpha",   # used in old BDT
    "mm_kin_alphaXY",
    "mm_kin_spvip",   # used in old BDT
    "mm_kin_pvip",   
    "mm_iso",         # used in old BDT
    ### good
    "mm_m1iso",       # used in old BDT
    "mm_m2iso",       # used in old BDT
    "mm_kin_sl3d",    # used in old BDT
    "mm_kin_vtx_chi2dof", # used in old BDT
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



if __name__ == "__main__":

    data = prepare_dataset()

    # x_train, x_test, y_train, y_test = get_train_and_test_datasets_split_randomly(data)
    event_index = 0
    x_train, x_test, y_train, y_test = get_train_and_test_datasets_split_be_event_number(data,3,event_index)

    print "Number of signal/background events in training sample: %u/%u (%0.3f)" % (sum(y_train==True),
                                                                                    sum(y_train==False),
                                                                                    sum(y_train==True)/float(sum(y_train==False)) )
    model_name = "Run2017-2018-%s-Event%u" % (datetime.now().strftime("%Y%m%d-%H%M"), event_index)
    
    # cross_validate(x_train, y_train, model_name)
    bst,dtrain,dtest = train(x_train, x_test, y_train, y_test,model_name)
    make_validation_plots(bst,dtrain,dtest,model_name)

