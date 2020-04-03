# import os
import sys
import pickle

import numpy as np
import uproot
from sklearn.model_selection import train_test_split


# from tqdm import tqdm
import xgboost as xgb

# from sklearn.metrics import roc_curve, roc_auc_score

# from utils import write_json, json_to_cfunc

from pprint import pprint

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

def _select_events(data,condition):
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
    train_data = _select_events(data,data['evt_event']%n!=index)
    test_data  = _select_events(data,data['evt_event']%n==index)
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
    param['eta'] = 0.3
    # max_depth [default=6]
    #
    # - Maximum depth of a tree. Increasing this value will make the model
    # more complex and more likely to overfit. 0 is only accepted in
    # lossguided growing policy when tree_method is set as hist and it
    # indicates no limit on depth. Beware that XGBoost aggressively
    # consumes memory when training a deep tree.  
    # - range: [0,Inf] (0 is only accepted in lossguided growing policy
    # when tree_method is set as hist)
    param['max_depth'] = 3
    param['silent'] = 1
    param['nthread'] = 15
    # eval_metric [default according to objective]
    # 
    # - Evaluation metrics for validation data, a default metric will be
    # assigned according to objective (rmse for regression, and error for
    # classification, mean average precision for ranking)
    # - auc: Area under the curve
    param['eval_metric'] = "auc"
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
    param['min_child_weight'] = 1.0
    # colsample_bytree, colsample_bylevel, colsample_bynode [default=1]
    # 
    # - This is a family of parameters for subsampling of columns.
    # - All colsample_by* parameters have a range of (0, 1], the default
    # value of 1, and specify the fraction of columns to be subsampled.
    # - colsample_bytree is the subsample ratio of columns when constructing
    # each tree. Subsampling occurs once for every tree constructed.
    param['colsample_bytree'] = 1.0
    
    return param

def train(x_train,x_test,y_train,y_test):
    """Train a model."""
    dtrain = xgb.DMatrix( x_train, label=y_train)
    dtest = xgb.DMatrix( x_test, label=y_test)
    evallist  = [(dtrain,'train'), (dtest,'eval')]

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
    param["scale_pos_weight"] = sumw_neg/float(sumw_pos)
    pprint(param)

    pklname = "train.pkl"
    bst = xgb.train( param.items(), dtrain, num_round, evallist, early_stopping_rounds=10 )
    bst.feature_names = feature_names
    max_length = max([len(s) for s in feature_names])
    feature_format = "%"+"%u"%max_length+"s"
    scores = bst.get_score(importance_type='gain')
    print "Importance scores:"
    for name in feature_names:
        if name in scores:
            print (feature_format+" %0.1f") % (name,scores[name])
        else:
            print (feature_format+" unused") % name
    pickle.dump(bst,open(pklname,"wb"))
    # write_json("model.json", bst, feature_names)
    # json_to_cfunc("model.json", fname_out="func.h")
    
feature_names = [
    #
    # Features in order of decreasing contribution
    #
    ### the best
    "mm_kin_alpha",   # used in old BDT
    "mm_kin_alphaXY",
    "mm_kin_spvip",   # used in old BDT
    "mm_iso",         # used in old BDT
    ### good
    "mm_m1iso",       # used in old BDT
    "mm_m2iso",       # used in old BDT
    "mm_kin_sl3d",    # used in old BDT
    ### weak 
    "mm_nBMTrks",    
    "mm_kin_vtx_chi2dof", # used in old BDT
    # "mm_mu2_pt",
    "mm_closetrks1",
    "mm_nDisTrks",
    ### very weak (in order of decreasing contribution)
    # "mm_mu2_eta",
    # "mm_docatrk",     # used in old BDT
    # "mm_closetrk",    # used in old BDT
    # "mm_kin_eta",     # used in old BDT
    ### new 
]

signal     = load_data("python/BsToMuMu.root")
background = load_data("python/QCD_Pt-80to120_MuEnrichedPt5_RunIIAutumn18MiniAOD.root")

data = merge_data([signal,background])

x_train, x_test, y_train, y_test = get_train_and_test_datasets_split_randomly(data)
# x_train, x_test, y_train, y_test = get_train_and_test_datasets_split_be_event_number(data)
train(x_train, x_test, y_train, y_test)

