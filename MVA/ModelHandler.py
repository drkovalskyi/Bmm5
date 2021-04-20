import sys, os, subprocess, json
from datetime import datetime

import numpy as np
import uproot
from sklearn.model_selection import train_test_split

import xgboost as xgb
from sklearn.metrics import roc_curve, roc_auc_score

from pprint import pprint

import matplotlib as mpl
# https://matplotlib.org/faq/usage_faq.html
mpl.use('Agg')
import matplotlib.pyplot as plt

class ModelHandler(object):
    """XGBoost Model training and validation"""
    def __init__(self, output_path, features, tree_name, event_branch_name):
        self.output_path = output_path
        self.features = features
        self.data = None
        self.tree_name = tree_name
        self.event_branch_name = event_branch_name
        self.x_train = None
        self.x_test  = None
        self.y_train = None
        self.y_test  = None
        self.bst = None
        self.train_matrix = None
        self.test_matrix = None

    def get_true_classification(self, data):
        """Make a column representing the true classification of events.

        The collumn is used for training and performance testing.
        """
        raise Exception("Not implemented")

    def load_data(self, file_name):
        """Load ROOT data and turn MVA tree into an array.

        The array is represented in a form of dictionary with the branch
        name representing a key and the value is numpy.ndarray.
        """
        print "Loading data from", file_name
        f = uproot.open(file_name)
        tree = f[self.tree_name]
        return tree.arrays(tree.keys())

    def merge_data(self, datasets):
        """Merge multiple datasets into one dataset."""
        dataset = dict()
        if len(datasets) < 1:
            return dataset
        branches = datasets[0].keys()
        for branch in branches:
            dataset[branch] = np.concatenate(
                [ds[branch] for ds in datasets],
                axis=None)
        return dataset

    def load_datasets(self, input_list):
        """Load and merge ROOT trees with MVA data into a single dataset."""
        datasets = []
        if self.data:
            datasets.append(self.data)
        for entry in input_list:
            files = subprocess.check_output("find %s -type f -name '*root'" % entry, shell=True)
            for f in files.splitlines():
                datasets.append(self.load_data(f))
        self.data = self.merge_data(datasets)
        print "Total number of events:", len(self.data[self.event_branch_name])

    def apply_selection(self, selection_mask):
        """Select events for the model using the selection mask"""
        self.data = self.select_events(self.data, selection_mask)
        print "Number of events selected:", len(self.data[self.event_branch_name])
        
    def get_feature_data(self, data):
        """Get columns containing feature data."""
        x_data = []
        for feature in self.features:
            x_data.append(data[feature])
        return np.column_stack(x_data)

    def select_events(self, data, condition):
        """Select events"""
        selected_data = dict()
        indices = np.where(condition)[0]
        for branch in data:
            selected_data[branch] = data[branch].take(indices)
        return selected_data

    def prepare_train_and_test_datasets(self, n=4, index=0):
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
        train_data = self.select_events(self.data, self.data[self.event_branch_name] % n != index)
        test_data  = self.select_events(self.data, self.data[self.event_branch_name] % n == index)
        self.y_train = self.get_true_classification(train_data)
        self.y_test  = self.get_true_classification(test_data)
        self.x_train = self.get_feature_data(train_data)
        self.x_test  = self.get_feature_data(test_data)


    def get_parameters(self):
        """Get parameters for XGBoost."""
        param = {}

        # https://xgboost.readthedocs.io/en/latest/parameter.html
        #
        # param['tree_method'] = 'hist'
        #

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


    def train(self, model_name="test", num_boost_round=5000, make_validation_plots=True):
        """Train a model"""
        self.train_matrix = xgb.DMatrix(self.x_train, label=self.y_train, feature_names=self.features)
        self.test_matrix  = xgb.DMatrix(self.x_test,  label=self.y_test,  feature_names=self.features)
        evallist  = [(self.train_matrix, 'train'), (self.test_matrix, 'eval')]
        param = self.get_parameters()

        sumw_pos = np.abs(self.train_matrix.get_label() == 1).sum()
        sumw_neg = np.abs(self.train_matrix.get_label() == 0).sum()
        # scale_pos_weight [default=1]
        #
        # - Control the balance of positive and negative weights, useful for
        # unbalanced classes. A typical value to consider: sum(negative
        # instances) / sum(positive instances). See Parameters Tuning for more
        # discussion.
        # https://xgboost.readthedocs.io/en/latest/tutorials/param_tuning.html

        # target_s_over_b = 0.005
        # target_s_over_b = .1

        # param["scale_pos_weight"] = sumw_neg/float(sumw_pos)*target_s_over_b
        
        pprint(param)

        self.bst = xgb.train(param.items(), self.train_matrix, num_boost_round, evallist, early_stopping_rounds=10)
        
        # bst.feature_names = self.features
        max_length = max([len(s) for s in self.features])
        feature_format = "%" + "%u" % max_length + "s"
        scores = self.bst.get_score(importance_type='gain')
        print "Importance scores:"
        for name in self.features:
            if name in scores:
                print (feature_format + " %0.1f") % (name, scores[name])
            else:
                print (feature_format + " unused") % name
        # pickle.dump({'model':bst,'feature_names':self.features},open("%s/%s.pkl"%(output_path,model_name),"wb"))
        if not os.path.exists(self.output_path):
            subprocess.call("mkdir -p %s" % self.output_path, shell=True)
        self.bst.dump_model("%s/%s.txt" % (self.output_path, model_name))
        self.bst.save_model("%s/%s.model" % (self.output_path, model_name))
        json.dump(self.features,
                  open("%s/%s.features" % (self.output_path, model_name), "w"))
        json.dump(param,
                  open("%s/%s.params" % (self.output_path, model_name), "w"))

        # write_json("model.json", bst, self.features)
        # json_to_cfunc("model.json", fname_out="func.h")
        if make_validation_plots:
            self.make_validation_plots(model_name)

    def make_validation_plots(self, model_name="test"):
        """Plot ROC curves and score distributions"""
        y_pred_train = self.bst.predict(self.train_matrix)
        y_pred_test  = self.bst.predict(self.test_matrix)

        fig, ax = plt.subplots()
        preds_bkg_test  = y_pred_test[self.y_test == 0]
        preds_sig_test  = y_pred_test[self.y_test == 1]
        preds_bkg_train = y_pred_train[self.y_train == 0]
        preds_sig_train = y_pred_train[self.y_train == 1]

        density = True
        bins = np.linspace(0.0, 1, 50)
        edges = 0.5 * (bins[:-1] + bins[1:])

        from scipy import stats
        pks_bkg = stats.ks_2samp(preds_bkg_train, preds_bkg_test)[1]
        pks_sig = stats.ks_2samp(preds_sig_train, preds_sig_test)[1]

        counts_train_bkg, _, _ = ax.hist(
            preds_bkg_train,
            bins=bins,
            histtype="stepfilled",
            alpha=0.45,
            normed=density,
            label="bkg, train",
            color=["b"])
        counts_train_sig, _, _ = ax.hist(
            preds_sig_train,
            bins=bins,
            histtype="stepfilled",
            alpha=0.45,
            normed=density,
            label="sig, train",
            color=["r"])
        counts_test_bkg, _, _ = ax.hist(
            preds_bkg_test,
            bins=bins,
            histtype="step",
            alpha=1.0,
            normed=density,
            label="bkg, test (KS prob = {:.2f})".format(pks_bkg),
            color=["b"],
            lw=1.5,
            linestyle="solid")
        counts_test_sig, _, _ = ax.hist(
            preds_sig_test,
            bins=bins,
            histtype="step",
            alpha=1.0,
            normed=density,
            label="sig, test (KS prob = {:.2f})".format(pks_sig),
            color=["r"],
            lw=1.5,
            linestyle="solid")

        ax.set_yscale("log")
        ax.legend()

        ax.set_ylim([0.01, ax.get_ylim()[1]])
        fig.set_tight_layout(True)
        fig.savefig("%s/%s-validation-disc.pdf" % (self.output_path, model_name))

        fpr_test, tpr_test, _ = roc_curve(self.y_test, y_pred_test)
        fpr_train, tpr_train, _ = roc_curve(self.y_train, y_pred_train)
        auc_test = roc_auc_score(self.y_test, y_pred_test)
        auc_train = roc_auc_score(self.y_train, y_pred_train)
        fig, ax = plt.subplots()
        plt.grid(which='both', axis='both')
        ax.plot(fpr_test, tpr_test, label="test AUC = {:.6f}".format(auc_test))
        ax.plot(fpr_train, tpr_train, label="train AUC = {:.6f}".format(auc_train))
        ax.set_xlabel("bkg eff")
        ax.set_ylabel("sig eff")
        ax.set_xlim([0.001, 1])
        ax.set_ylim([0.2, 1.0])
        ax.set_title("ROC curves")
        ax.legend()
        fig.set_tight_layout(True)
        ax.set_yscale("log")
        ax.set_xscale("log")
        fig.savefig("%s/%s-validation-roc.pdf" % (self.output_path, model_name))




if __name__ == "__main__":

    feature_names = [
        ### Good
        "trkValidFrac",
        "glbTrackProbability",
        "nLostHitsInner",
        "nLostHitsOuter",
        "trkKink",
        "chi2LocalPosition",
        "match2_dX",
        "match2_pullX",
        "match1_dX",
        "match1_pullX",

        ### Weak but useful
        "nPixels",
        "nValidHits",
        "nLostHitsOn",
        "eta",
        "match2_dY",
        "match1_dY",
        "match2_pullY",
        "match1_pullY",
        "match2_pullDyDz",
        "match1_pullDyDz",
        "match2_pullDxDz",
        "match1_pullDxDz",

    ]

    files = [
        # "/eos/cms/store/group/phys_bphys/bmm/bmm5/PostProcessing/FlatNtuples/512/muon_mva/QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_13TeV_pythia8+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v3+MINIAODSIM/15475f7d44f873b2f0aaca4372b12284.root"
        "/eos/cms/store/group/phys_bphys/bmm/bmm5/PostProcessing/FlatNtuples/512/muon_mva/QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_13TeV_pythia8+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v3+MINIAODSIM/"
    ]

    class MuonMVA(ModelHandler):
        def get_true_classification(self, data):
            return (data["sim_type"] != 1)
        def get_parameters(self):
            parameters = super(MuonMVA, self).get_parameters()
            parameters['min_child_weight'] = 0.01
            parameters['max_depth'] = 5
            parameters['eta'] = 1
            parameters['tree_method'] = 'hist' 
            return parameters
        
    model = MuonMVA("results/test", feature_names, "muons", "evt")
    model.load_datasets(files)
    model.apply_selection(model.data["pt"]>5)

    event_index = 0

    model.prepare_train_and_test_datasets(3, event_index)
    
    print "Number of signal/background events in training sample: %u/%u (%0.3f)" \
        % (sum(model.y_train == True),
           sum(model.y_train == False),
           sum(model.y_train == True) / float(sum(model.y_train == False)))

    model_name = "Run2018-%s-Event%u" % (datetime.now().strftime("%Y%m%d-%H%M"), event_index)

    # model.cross_validate(model_name)
    model.train(model_name)
