import ROOT
from sklearn.metrics import roc_curve, roc_auc_score
import os
import matplotlib
matplotlib.use('Agg')  # Use a non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns
# Set the color palette
sns.set_palette("colorblind")  # Or "Set2", "tab10", etc.
from scipy.interpolate import interp1d

limit = 10000
output_path = "/afs/cern.ch/user/d/dmytro/www/public_html/plots/muonid_run3"
nanoaod_path = "/eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/527/"
skim_path = "/eos/cms/store/group/phys_bphys/bmm/bmm6/PostProcessing/Skims/527/tau3mu"
min_pt = 2.0
max_pt = 9999

samples = {
    'tau3mu':[
	nanoaod_path + "/DstoTau_Tauto3Mu_3MuFilter_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v2+MINIAODSIM/*root",
    ],
    'tau3mu_skim':[
	skim_path + "/DstoTau_Tauto3Mu_3MuFilter_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v2+MINIAODSIM/*root",
    ],
    'inclusive':[
        nanoaod_path + "/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8+Run3Summer22MiniAODv3-Pilot_124X_mcRun3_2022_realistic_v12-v5+MINIAODSIM/1a*root",
    ],
    'inclusive_skim':[
        skim_path + "/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8+Run3Summer22MiniAODv3-Pilot_124X_mcRun3_2022_realistic_v12-v5+MINIAODSIM/*root",
    ],
    'data_skim':[
        skim_path + "/ParkingDoubleMuonLowMass0+Run2022C-PromptReco-v1+MINIAOD/*root",
    ],
}

regions = {
    # barrel
    # 'eta_0.0-0.8_pt_all': {
    #     'eta_min':0.0,
    #     'eta_max':0.8,
    #     'pt_min':0,
    #     'pt_max':9999,
    #     'title': '$|\eta|<0.8$'
    # },
    # 'eta_0.0-0.8_pt_0-4': {
    #     'eta_min':0.0,
    #     'eta_max':0.8,
    #     'pt_min':0,
    #     'pt_max':4,
    #     'title': '$|\eta|<0.8$ and $p_T<4$'
    # },
    'eta_0.0-0.8_pt_4-6': {
        'eta_min':0.0,
        'eta_max':0.8,
        'pt_min':4,
        'pt_max':6,
        'title': '$|\eta|<0.8$ and $p_T\in[4,6]$'
    },
    'eta_0.0-0.8_pt_6-10': {
        'eta_min':0.0,
        'eta_max':0.8,
        'pt_min':6,
        'pt_max':10,
        'title': '$|\eta|<0.8$ and $p_T\in[6,10]$'
    },
    # 'eta_0.0-0.8_pt_10-inf': {
    #     'eta_min':0.0,
    #     'eta_max':0.8,
    #     'pt_min':10,
    #     'pt_max':9999,
    #     'title': '$|\eta|<0.8$ and $p_T>10$'
    # },
    
    # # overlap region
    # 'eta_0.8-1.2_pt_all': {
    #     'eta_min':0.8,
    #     'eta_max':1.2,
    #     'pt_min':0,
    #     'pt_max':9999,
    #     'title': '$|\eta|\in[0.8,1.2]$'
    # },
    # 'eta_0.8-1.2_pt_0-4': {
    #     'eta_min':0.8,
    #     'eta_max':1.2,
    #     'pt_min':0,
    #     'pt_max':4,
    #     'title': '$|\eta|\in[0.8,1.2]$ and $p_T<4$'
    # },
    'eta_0.8-1.2_pt_4-6': {
        'eta_min':0.8,
        'eta_max':1.2,
        'pt_min':4,
        'pt_max':6,
        'title': '$|\eta|\in[0.8,1.2]$ and $p_T\in[4,6]$'
    },
    'eta_0.8-1.2_pt_6-10': {
        'eta_min':0.8,
        'eta_max':1.2,
        'pt_min':6,
        'pt_max':10,
        'title': '$|\eta|\in[0.8,1.2]$ and $p_T\in[6,10]$'
    },
    # 'eta_0.8-1.2_pt_10-inf': {
    #     'eta_min':0.8,
    #     'eta_max':1.2,
    #     'pt_min':10,
    #     'pt_max':9999,
    #     'title': '$|\eta|\in[0.8,1.2]$ and $p_T>10$'
    # },
    
    # # endcap
    # 'eta_1.2-2.1_pt_all': {
    #     'eta_min':1.2,
    #     'eta_max':2.1,
    #     'pt_min':0,
    #     'pt_max':9999,
    #     'title': '$|\eta|\in[1.2,2.1]$'
    # },
    # 'eta_1.2-2.1_pt_0-4': {
    #     'eta_min':1.2,
    #     'eta_max':2.1,
    #     'pt_min':0,
    #     'pt_max':4,
    #     'title': '$|\eta|\in[1.2,2.1]$ and $p_T<4$'
    # },
    'eta_1.2-2.1_pt_4-6': {
        'eta_min':1.2,
        'eta_max':2.1,
        'pt_min':4,
        'pt_max':6,
        'title': '$|\eta|\in[1.2,2.1]$ and $p_T\in[4,6]$'
    },
    'eta_1.2-2.1_pt_6-10': {
        'eta_min':1.2,
        'eta_max':2.1,
        'pt_min':6,
        'pt_max':10,
        'title': '$|\eta|\in[1.2,2.1]$ and $p_T\in[6,10]$'
    },
    # 'eta_1.2-2.1_pt_10-inf': {
    #     'eta_min':1.2,
    #     'eta_max':2.1,
    #     'pt_min':10,
    #     'pt_max':9999,
    #     'title': '$|\eta|\in[1.2,2.1]$ and $p_T>10$'
    # },
}

LorentzVector = ROOT.ROOT.Math.LorentzVector('ROOT::Math::PtEtaPhiM4D<double>')
muon_mass = 0.10566

chains = dict()

def load_data(sample_name):
    chain = ROOT.TChain("Events")
    for sample in samples[sample_name]:
        chain.Add(sample)
            
    print("Total number of events in sample %s : %u" % (sample_name, chain.GetEntries()))
    return chain


def get_data(sample_name):
    if sample_name not in chains:
        chains[sample_name] = load_data(sample_name)
    return chains[sample_name]


def good_candidate(event, i):
    if not event.Muon_isGlobal[i]: return False
    if not event.Muon_isPFcand[i]: return False
    if event.Muon_pt[i] < min_pt: return False
    if event.Muon_pt[i] > max_pt: return False
    return True


def select_muons(chain, ids, type='all', region=None, trigger=False, trigger_cut=True,
                 preselection=False):
    """Select muons for muon id studies

    Parameters:

    * type: define if either signal ('signal') or background ('bkg')
      muon candidates to keep or both ('all')
    * region: dictionaty with pt and eta cuts
    * trigger: use trigger requirement
    * trigger_cut: apply trigger as a pre-selection cut or create a
      new array and set prediction probability to zero
    * preselection: apply muon preselection
    """
    results = {
        'label':[],
    }
    for id in ids:
        results[id] = []
        if trigger and not trigger_cut:
            results[id + "_trig"] = []

    n_sig_before_trigger = 0
    n_bkg_before_trigger = 0
    for event in chain:
        for i in range(event.nMuon):
            if preselection and not good_candidate(event,i):
                continue
            if region != None:
                cuts = regions[region]
                if abs(event.Muon_eta[i]) < cuts['eta_min'] or abs(event.Muon_eta[i]) > cuts['eta_max']:
                    continue
                if event.Muon_pt[i] < cuts['pt_min'] or event.Muon_pt[i] > cuts['pt_max']:
                    continue
                
            label = None
            if event.MuonId_simType[i] == 3:
                label = 1
            if event.MuonId_simType[i] == 2:
                label = 0
                
            if label == None or \
               (type == 'signal' and label == 0) or \
               (type == 'bkg' and label == 1):
                continue

            if label:
                n_sig_before_trigger += 1
            else:
                n_bkg_before_trigger += 1
                
            if trigger and trigger_cut:
                if event.MuonId_hlt_pt[i] < 0:
                    continue

            results['label'].append(label)
            for id in ids:
                results[id].append(getattr(event,id)[i])
                if trigger and not trigger_cut:
                    results[id + "_trig"].append(getattr(event,id)[i] if event.MuonId_hlt_pt[i] > 0 else 0)
                
        if limit and len(results['label']) >= limit:
            break
                
    print("Number of positive labels: ", sum(results['label']))
    print("Number of negative labels: ", len(results['label']) - sum(results['label']))
    if trigger:
        print("Number of positive labels before trigger: ", n_sig_before_trigger)
        print("Number of negative labels before trigger: ", n_bkg_before_trigger)
        
    return results


def select_muons_for_tau3mu(chain, ids, type='all', preselection=True, dimuon_trigger=False):
    """Select muons for muon id studies with tau3mu analysis selections

    Parameters:

    * type: define if either signal ('signal') or background ('bkg')
      muon candidates to keep or both ('all'). Also support labeling
      all events either signal or background with any type matching:
      'label_signal' and 'label_bkg'. In addition to that support
      'label_sig_as_bkg'.
    * preselection: apply muon preselection
    
    * dimuon_trigger: apply HLT_DoubleMu4_3_LowMass trigger to the
      event and require the first two muons to be triggering.

    """
    if type not in ['all', 'signal', 'bkg', 'label_signal', 'label_bkg',
                    'label_sig_as_bkg']:
        raise Exception("Unknown type: %s" % type)
    
    results = {
        'label':[],
    }
    for id in ids:
        results[id] = []

    sim_types = dict()
    for event in chain:
        if event.nMuon < 3:
            continue

        if dimuon_trigger and not event.HLT_DoubleMu4_3_LowMass:
            continue
        
        muon = None
        for mu1 in range(event.nMuon - 2):
            if preselection and not good_candidate(event, mu1):
                continue
            if dimuon_trigger and event.MuonId_hlt_pt[mu1] < 0:
                continue
            p1 = LorentzVector(event.Muon_pt[mu1], event.Muon_eta[mu1], event.Muon_phi[mu1], muon_mass)
            for mu2 in range(mu1 + 1, event.nMuon - 1):
                if preselection and not good_candidate(event, mu2):
                    continue
                if dimuon_trigger and event.MuonId_hlt_pt[mu2] < 0:
                    continue
                p2 = LorentzVector(event.Muon_pt[mu2], event.Muon_eta[mu2], event.Muon_phi[mu2], muon_mass)
                for mu3 in range(mu2 + 1, event.nMuon):
                    if preselection and not good_candidate(event, mu3):
                        continue
                    p3 = LorentzVector(event.Muon_pt[mu3], event.Muon_eta[mu3], event.Muon_phi[mu3], muon_mass)
                    tau_mass = (p1 + p2 + p3).mass()
                    if tau_mass > 2.0 or tau_mass < 1.6:
                        continue
                    muon = mu3
                    break
                if muon != None:
                    break
            if muon != None:
                break

        if muon == None:
            continue

        if hasattr(event, "MuonId_simType"):
            if event.MuonId_simType[muon] not in sim_types:
                sim_types[event.MuonId_simType[muon]] = 0
            sim_types[event.MuonId_simType[muon]] += 1
        
        label = None

        if type == 'label_signal':
            label = 1
        elif type == 'label_bkg':
            label = 0
        else:
            if event.MuonId_simType[muon] == 3:
                if type == 'all' or type == 'signal':
                    label = 1
                elif type == 'label_sig_as_bkg':
                    label = 0
            if event.MuonId_simType[muon] == 2:
                if type == 'all' or type == 'bkg':
                    label = 0
                
        if label == None:
            continue

        results['label'].append(label)
        for id in ids:
            results[id].append(getattr(event,id)[muon])
                
        if limit and len(results['label']) >= limit:
            break

    print("Number of positive labels: ", sum(results['label']))
    print("Number of negative labels: ", len(results['label']) - sum(results['label']))

    print("3d muon by sim-type distribution:")
    for sim_type, count in sorted(sim_types.items(), key=lambda x: x[1], reverse=True):
        print("\t%d: \t%u" % (sim_type, count))
    
    return results


def select_muons_for_pseudo_fake_study(chain, ids, preselection=False, min_dpt=0):
    """Pseudo background study

    Select heavy flavor muon pairs and label higher pt one as signal
    and the low pt one as background

    """
    results = {
        'label':[],
    }
    for id in ids:
        results[id] = []
        
    for event in chain:
        if event.nMuon < 2:
            continue

        muon1 = None
        muon2 = None
        for mu1 in range(event.nMuon - 1):
            if preselection and not good_candidate(event, mu1):
                continue
            if not event.MuonId_simType[mu1] == 3:
                continue
            for mu2 in range(mu1 + 1, event.nMuon):
                if preselection and not good_candidate(event, mu2):
                    continue
                if not event.MuonId_simType[mu2] == 3:
                    continue
                if abs(event.Muon_pt[mu1] - event.Muon_pt[mu2]) < min_dpt:
                    continue
                if event.Muon_pt[mu1] > event.Muon_pt[mu2]:
                    muon1 = mu1
                    muon2 = mu2
                else:
                    muon1 = mu2
                    muon2 = mu1
                break

        if muon1 != None and muon2 != None:

            results['label'].append(1)
            for id in ids:
                results[id].append(getattr(event,id)[muon1])
                
            results['label'].append(0)
            for id in ids:
                results[id].append(getattr(event,id)[muon2])
                
        if limit and len(results['label']) >= limit:
            break

    print("Number of positive labels: ", sum(results['label']))
    print("Number of negative labels: ", len(results['label']) - sum(results['label']))
    return results


def plot_roc_curve(data, ids, filename,
                   xlabel="Background Efficiency", ylabel="Signal Efficiency",
                   title="Receiver Operating Characteristic",
                   test_fpr = None):
    plt.figure(figsize=(6,6), dpi=150)

    points = []
    for id in ids:
        fpr, tpr, thresholds = roc_curve(data['label'], data[id])
        roc_line, = plt.plot(fpr, tpr, lw=2, label=id)

        if test_fpr != None:
            # Interpolate
            interp_tpr = interp1d(fpr, tpr, kind='nearest')
            estimated_tpr = interp_tpr(test_fpr)

            # Extract the color of the line
            line_color = roc_line.get_color()

            points.append((test_fpr, estimated_tpr, line_color))

    if len(points) > 0:
        points.sort(key=lambda x: x[1], reverse=True)
        offset_x = 0.02
        for i, (test_fpr, estimated_tpr, line_color) in enumerate(points):
            plt.scatter(test_fpr, estimated_tpr, color=line_color)
            if i == 0:
                offset_y = 0.05
                offset_x = -0.1
            elif i == len(points) - 1:
                offset_y = -0.05
                offset_x = 0.02
            else:
                offset_y = 0
                offset_x = 0.02
            plt.text(test_fpr + offset_x, estimated_tpr + offset_y, f'({estimated_tpr:.3f})', color=line_color)
    
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend(loc="lower right")

    if not os.path.exists(output_path):
        os.makedirs(output_path)

    plt.savefig("%s/%s.png" % (output_path, filename))
    plt.savefig("%s/%s.pdf" % (output_path, filename))

def merge_data(data1, data2):
    result = {}
    for entry in data1:
        result[entry] = data1[entry] + data2[entry]
    return result
    

ids = [
    "MuonId_xgb_pt2",
    "MuonId_xgb_weighted_pt2",
    "MuonId_ewkMvaId",
    "Muon_softMva"
]

# # signal = select_muons(load_data('tau3mu'), ids, "signal")
# # bkg = select_muons(load_data('tau3mu'), ids, "bkg")
# # inclusive_bkg = select_muons(load_data('inclusive'), ids, "bkg")

# # signal = select_muons_for_tau3mu(load_data('tau3mu_skim'), ids, 1)
# signal = select_muons_for_tau3mu(load_data('tau3mu_skim'), ids)
# print(len(signal['label']))

# # bkg = select_muons_for_tau3mu(load_data('inclusive_skim'), ids, 0)
# bkg = select_muons_for_tau3mu(load_data('inclusive_skim'), ids)
# print(len(bkg['label']))

## All IDs
# data = select_muons(get_data('inclusive'), ids, "all")
#
# plot_roc_curve(data, ids, "inclusive_reference_all",
#                "Decay in flight efficiency", "Heavy flavour muon efficiency",
#                "Inclusive Dilepton MinBias Events")

ids = [
    # "MuonId_xgb_pt2",
    "MuonId_xgb_weighted_pt2",
    "MuonId_ewkMvaId",
    "Muon_softMva"
]

## Relevant IDs
# plot_roc_curve(data, ids, "inclusive_reference",
#                "Decay in flight efficiency", "Heavy flavour muon efficiency",
#                "Inclusive Dilepton MinBias Events")

# ## Performace in different regions
# for region in regions:
#     data = select_muons(get_data('inclusive'), ids, "all", region)

#     plot_roc_curve(data, ids, "inclusive_reference" + region,
#                    "Decay in flight efficiency", "Heavy flavour muon efficiency",
#                    "Inclusive Dilepton MinBias Events (%s)" % regions[region]['title'])

## Trigger
# data = select_muons(get_data('inclusive'), ids, "all", None, True)

# plot_roc_curve(data, ids, "inclusive_triggered",
#                "Decay in flight efficiency", "Heavy flavour muon efficiency",
#                "Inclusive Dilepton MinBias Events (triggered)",
#                0.5)

# # Trigger overall impact
# ids = ["MuonId_xgb_weighted_pt2"]
# data = select_muons(get_data('inclusive'), ids, "all", None,
#                     trigger=True, trigger_cut=False)

# ids = [
#     "MuonId_xgb_weighted_pt2",
#     "MuonId_xgb_weighted_pt2_trig",
# ]
# plot_roc_curve(data, ids, "inclusive_trigger_impact",
#                "Decay in flight efficiency", "Heavy flavour muon efficiency",
#                "Inclusive Dilepton MinBias Events",
#                0.3)

#############################################
#         Triger regional study
#############################################
# for region in regions:
#     ids = ["MuonId_xgb_weighted_pt2"]
#     data = select_muons(get_data('inclusive'), ids, "all", region,
#                         trigger=True, trigger_cut=False)

#     ids = [
#         "MuonId_xgb_weighted_pt2",
#         "MuonId_xgb_weighted_pt2_trig",
#     ]
#     plot_roc_curve(data, ids, "inclusive_trigger_impact" + region,
#                    "Decay in flight efficiency", "Heavy flavour muon efficiency",
#                    "Inclusive Dilepton MinBias Events (%s)" % regions[region]['title'],
#                    0.3)
    

#############################################
#         tau3mu studies
#############################################
# ## muon object level comparison
# signal = select_muons(get_data('tau3mu'), ids, "signal")
# bkg = select_muons(get_data('inclusive'), ids, "bkg")
# data = merge_data(signal, bkg)

# ids = [
#     "MuonId_xgb_weighted_pt2",
#     "MuonId_ewkMvaId",
#     "Muon_softMva"
# ]
# plot_roc_curve(data, ids, "tau3mu_vs_inclusive-just_muons",
#                "Decay in flight efficiency", "Heavy flavour muon efficiency",
#                r'Signal $\tau\to3\mu$ vs Background Inclusive Dilepton MinBias')

# # tau3mu preselection impact
# data = select_muons(get_data('inclusive'), ids, preselection=True)
# plot_roc_curve(data, ids, "inclusive_tau3mu_preselection",
#                "Decay in flight efficiency", "Heavy flavour muon efficiency",
#                "Inclusive Dilepton MinBias Events (tau3mu peselection)",
#                0.3)

# # tau3mu preselection+trigger impact
# data = select_muons(get_data('inclusive'), ids, preselection=True,
#                     trigger=True)
# plot_roc_curve(data, ids, "inclusive_tau3mu_preselection_n_trigger",
#                "Decay in flight efficiency", "Heavy flavour muon efficiency",
#                r"Inclusive Dilepton MinBias Events ($\tau\to3\mu$ peselection + trigger)",
#                0.3)

# # tau3mu preselection+trigger impact by region
# for region in regions:
#     data = select_muons(get_data('inclusive'), ids, preselection=True,
#                         region=region, trigger=True)
#     plot_roc_curve(data, ids, "inclusive_tau3mu_preselection_n_trigger_" + region,
#                    "Decay in flight efficiency", "Heavy flavour muon efficiency",
#                    r"Inclusive ($\tau\to3\mu$ peselection + trig + %s)" % regions[region]['title'])

# # # psuedo fake study
# data = select_muons_for_pseudo_fake_study(get_data('tau3mu'), ids, min_dpt=1)
# plot_roc_curve(data, ids, "tau3mu_pseudo_fakes_min_dpt_1",
#                "Low $p_T$ muon efficiency", "Hight $p_T$ muon efficiency",
#                r"Signal muons from $\tau\to3\mu$")

## 3d muon - not matched
# signal = select_muons_for_tau3mu(get_data('tau3mu_skim'), ids, type='label_signal', preselection=False)
# bkg = select_muons_for_tau3mu(get_data('inclusive_skim'), ids, type='label_bkg', preselection=False)
# data = merge_data(signal, bkg)
# plot_roc_curve(data, ids, "tau3mu_vs_inclusive-3d_muon-no_match-no_preselection",
#                "Inclusive dilepton MinBias efficiency", r"$\tau\to3\mu$ efficiency",
#                r'3d muon efficiency (no type matching, no preselection)')
# signal = select_muons_for_tau3mu(get_data('tau3mu_skim'), ids, type='label_signal', preselection=True)
# bkg = select_muons_for_tau3mu(get_data('inclusive_skim'), ids, type='label_bkg', preselection=True)
# data = merge_data(signal, bkg)
# plot_roc_curve(data, ids, "tau3mu_vs_inclusive-3d_muon-no_match-preselection",
#                "Inclusive dilepton MinBias efficiency", r"$\tau\to3\mu$ efficiency",
#                r'3d muon efficiency (no type matching, preselection)')
# Trigger
# signal = select_muons_for_tau3mu(get_data('tau3mu_skim'), ids, type='label_signal',
#                                  preselection=False, dimuon_trigger=True)
# bkg = select_muons_for_tau3mu(get_data('inclusive_skim'), ids, type='label_bkg',
#                               preselection=False, dimuon_trigger=True)
# data = merge_data(signal, bkg)
# plot_roc_curve(data, ids, "tau3mu_vs_inclusive-3d_muon-no_match-no_preselection-dimuon_trig",
#                "Inclusive dilepton MinBias efficiency", r"$\tau\to3\mu$ efficiency",
#                r'3d muon efficiency (no matching, no preselection, dimuon trigger)')
# signal = select_muons_for_tau3mu(get_data('tau3mu_skim'), ids, type='label_signal',
#                                  preselection=True, dimuon_trigger=True)
# bkg = select_muons_for_tau3mu(get_data('inclusive_skim'), ids, type='label_bkg',
#                               preselection=True, dimuon_trigger=True)
#
# data = merge_data(signal, bkg)
# plot_roc_curve(data, ids, "tau3mu_vs_inclusive-3d_muon-no_match-preselection-dimuon_trig",
#                "Inclusive dilepton MinBias efficiency", r"$\tau\to3\mu$ efficiency",
#                r'3d muon efficiency (no matching, preselection, dimuon trigger)')

# # Real data
signal = select_muons_for_tau3mu(get_data('tau3mu_skim'), ids, type='label_signal',
                                 preselection=True, dimuon_trigger=True)
bkg = select_muons_for_tau3mu(get_data('data_skim'), ids, type='label_bkg',
                              preselection=True, dimuon_trigger=True)
data = merge_data(signal, bkg)
plot_roc_curve(data, ids, "tau3mu_vs_data-3d_muon-no_match-preselection-dimuon_trig",
               "Run2022C efficiency", r"$\tau\to3\mu$ efficiency",
               r'3d muon efficiency (preselection, dimuon trigger)')

# Real data
# signal = select_muons_for_tau3mu(get_data('tau3mu_skim'), ids, type='label_signal',
#                                  preselection=False, dimuon_trigger=True)
# bkg = select_muons_for_tau3mu(get_data('data_skim'), ids, type='label_bkg',
#                               preselection=False, dimuon_trigger=True)
# data = merge_data(signal, bkg)
# plot_roc_curve(data, ids, "tau3mu_vs_data-3d_muon-no_preselection-dimuon_trig",
#                "Run2022C efficiency", r"$\tau\to3\mu$ efficiency",
#                r'3d muon efficiency (no preselection, dimuon trigger)')


# # 3d muon - signal-signal
# signal = select_muons_for_tau3mu(get_data('tau3mu_skim'), ids, type='signal', preselection=False)
# bkg = select_muons_for_tau3mu(get_data('inclusive_skim'), ids, type='label_sig_as_bkg', preselection=False)
# data = merge_data(signal, bkg)
# plot_roc_curve(data, ids, "tau3mu_vs_inclusive-3d_muon-sig-sig_no-preselection",
#                "Inclusive dilepton MinBias efficiency", r"$\tau\to3\mu$ efficiency",
#                r'3d muon efficiency (heavy flavor muons, no preselection)')


# signal = select_muons_for_tau3mu(get_data('tau3mu_skim'), ids, type='label_signal', preselection=False)
# bkg = select_muons_for_tau3mu(get_data('inclusive_skim'), ids, type='label_bkg', preselection=False)
# data = merge_data(signal, bkg)
# plot_roc_curve(data, ids, "tau3mu_vs_inclusive-3d_muon-no_match-no_preselection",
#                "Decay in flight efficiency", "Heavy flavour muon efficiency",
#                r'$\tau\to3\mu$ vs MinBias (3d muon, no matching, no preselection)')

