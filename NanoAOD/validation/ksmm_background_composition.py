import ROOT
import glob
import os

# input_path = '/data/dmytro/Run3-Bmm-NanoAODv12/root-files/InclusiveDileptonMinBias.root'
# input_path = '/tmp/dmytro/InclusiveDileptonMinBias.root'
# input_path = '/eos/cms/store/group/phys_bphys/bmm/bmm6/PostProcessing/Skims/528/ksmm/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8+Run3Summer22MiniAODv3-Pilot_124X_mcRun3_2022_realistic_v12-v5+MINIAODSIM/49d2e6c8092880787d26f917d0ba7fc2.root'
input_path = '/eos/cms/store/group/phys_bphys/bmm/bmm6/PostProcessing/Skims/530/ksmm/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8+Run3Summer22MiniAODv4-validDigi_130X_mcRun3_2022_realistic_v5-v4+MINIAODSIM/*.root'
# input_path = '/eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/526/K0sToMuMu_K0sFilter_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v1+MINIAODSIM/*.root'
# input_path = '/data/dmytro/Run3-Bmm-NanoAODv12/root-files/KsToMuMu.root'

verbose = True
# verbose_limit = 100
verbose_limit = None
origin_analysis = False
force_recreate = False
print_only_unclassified = True

pdgIds = {
    511: 'B0', 521: 'B+', 531: 'Bs', 541: 'Bc',
    411: 'D0', 421: 'D+', 431: 'Ds',
    321: 'K+', 211: 'pi+', 15: 'tau',
}

unknown_pdgIds = []


def name(pdgId):
    id = abs(pdgId)
    if id not in pdgIds:
        if id not in unknown_pdgIds:
            unknown_pdgIds.append(id)
        return "other"
    return pdgIds[id]


def print_decay(event, gen_index, indent):
    print("%spdgId: %d \tpt: %0.1f \teta: %0.2f \tstatus: %d" % (
        indent,
        event.GenPart_pdgId[gen_index],
        event.GenPart_pt[gen_index],
        event.GenPart_eta[gen_index],
        event.GenPart_status[gen_index]
        )
    )
    if event.GenPart_status[gen_index] != 1:
        for i in range(event.nGenPart):
            if event.GenPart_genPartIdxMother[i] == gen_index:
                print_decay(event, i, indent + "\t")

def hadron_type(id):

    digits = [int(digit) for digit in str(abs(id))]
    
    # Single b-quark
    # mesons
    if len(digits) >= 3 and digits[-3] == 5 and \
       (len(digits) == 3 or digits[-4] == 0):
        return "b-hadron"

    # Single c-quark
    # mesons
    if len(digits) >= 3 and digits[-3] == 4 and \
       (len(digits) == 3 or digits[-4] == 0):
        return "c-hadron"

    return None

def is_DtoMuNuX_decay(event, p_id, p_index):
    if hadron_type(p_id) != 'c-hadron':
        return False
    mu = 0
    nu = 0
    for gen_index in range(event.nGenPart):
        if event.GenPart_genPartIdxMother[gen_index] == p_index:
            id = event.GenPart_pdgId[gen_index]
            if hadron_type(id) == 'c-hadron':
                return is_DtoMuNuX_decay(event, id, gen_index)
            if abs(id) == 13:
                mu = id
            if abs(id) == 14:
                nu = id
    if mu * nu == -13 * 14:
        return True
    return False

def has_fake_muon(event, p_id, p_index):
    mu = 0
    nu = 0
    for gen_index in range(event.nGenPart):
        if event.GenPart_genPartIdxMother[gen_index] == p_index:
            id = event.GenPart_pdgId[gen_index]
            if abs(id) not in [13, 14]:
                if has_fake_muon(event, id, gen_index):
                    return True
            if abs(id) == 13:
                mu = id
            if abs(id) == 14:
                nu = id
    if mu * nu == -13 * 14 and abs(p_id) in [211, 321]:
        return True
    return False

def classify_event(event, ancestor, ancestor_index):

    # No common ancestor
    if ancestor == 0:
        return "No common ancestor"

    ancestor_type = hadron_type(ancestor)
    
    # Ks to pipi
    if abs(ancestor) == 310:
        signature = 1
        for gen_index in range(event.nGenPart):
            if event.GenPart_genPartIdxMother[gen_index] == ancestor_index:
                signature *= event.GenPart_pdgId[gen_index]
        if signature == -211 * 211:
            return "KsToPiPi"

    # B decay
    if ancestor_type == "b-hadron":
        mu = 0
        nu = 0
        n_charm_semileptonic = 0
        n_fakes = 0
        for gen_index in range(event.nGenPart):
            if event.GenPart_genPartIdxMother[gen_index] == ancestor_index:
                id = event.GenPart_pdgId[gen_index]
                if abs(id) == 13:
                    mu = id
                if abs(id) == 14:
                    nu = id
                if is_DtoMuNuX_decay(event, id, gen_index):
                    n_charm_semileptonic += 1
                if has_fake_muon(event, id, gen_index):
                    n_fakes += 1
        if mu * nu == -13 * 14:
            if n_charm_semileptonic == 1:
                return "H_b to H_c mu nu, H_c to mu nu X"
            if n_fakes == 1:
                return "H_b to X mu nu, X contains decay in flight"
        elif n_charm_semileptonic == 1:
            return "H_b to H_c X, H_c to mu nu X'"
        elif n_charm_semileptonic == 2:
            return "H_b to H_c H_c', H_c to mu nu X"

    # if ancestor_type == "b-hadron":
    #     mu = 0
    #     nu = 0
    #     chadron = 0
    #     for gen_index in range(event.nGenPart):
    #         if event.GenPart_genPartIdxMother[gen_index] == ancestor_index:
    #             id = event.GenPart_pdgId[gen_index]
    #             if abs(id) == 13:
    #                 mu = id
    #             if abs(id) == 14:
    #                 nu = id
    #             if hadron_type(id) == "c-hadron":
    #                 chadron = id
    #     if mu * nu == -13 * 14 and chadron != 0:
    #         return "H_b to H_c mu nu, H_c to mu nu X"

    # D semi-leptonic
    if is_DtoMuNuX_decay(event, ancestor, ancestor_index):
        return "H_c to mu nu X"

    return "Unclassified"


skim_file_name = "/tmp/dmytro/ksmm_skim.root"

if not os.path.exists(skim_file_name) or force_recreate:
    # Read data
    chain = ROOT.TChain("Events")
    # for f in glob.glob("%s/*.root" % input_path):
    #     chain.Add(f)
    chain.Add(input_path)

    # Skimming data
    df = ROOT.RDataFrame(chain)
    n_events = df.Count().GetValue()
    print("Number of events to process: %d" % n_events)

    df2 = df.Define("goodCandidates",
                    "mm_kin_mass>0.40 && mm_kin_mass<0.60 && mm_kin_sl3d > 3 && mm_kin_alpha < 0.1")
    df_final = df2.Filter("Sum(goodCandidates) > 0", "Event has good candidates")

    df_final.Snapshot("Events", skim_file_name, "nMuonId|MuonId_.*|nMuon|Muon_.*|nmm|mm_.*|nGenPart|GenPart_.*|event|run")

signatures = dict()

n_fakes = dict()

f = ROOT.TFile(skim_file_name)
events = f.Get("Events")

if not events:
    raise Exception("Cannot read the skim")

print("Number of events in the skim: %d" % events.GetEntries())

n_ref = 0
n_tight = 0
n_very_tight = 0
n_extreme = 0
verbose_count = 0

processes = dict()

# for event in chain:
for event in events:
    for mm_index in range(event.nmm):
        ## Selection
        if event.mm_kin_mass[mm_index] < 0.48 or event.mm_kin_mass[mm_index] > 0.60: continue
        if event.Muon_charge[event.mm_mu1_index[mm_index]] * \
           event.Muon_charge[event.mm_mu2_index[mm_index]] > 0: continue
        # if event.Muon_softMva[event.mm_mu1_index[mm_index]] < 0.45: continue
        # if event.Muon_softMva[event.mm_mu2_index[mm_index]] < 0.45: continue
        if event.mm_mu1_pt[mm_index] < 4: continue
        if event.mm_mu2_pt[mm_index] < 3: continue

        # reference
        if event.mm_kin_slxy[mm_index] < 10: continue
        if event.mm_kin_lxy[mm_index] < 1: continue
        if event.mm_kin_vtx_prob[mm_index] < 0.01: continue
        if event.mm_kin_alpha[mm_index] > 0.001: continue

        n_ref += 1
        
        # extra
        
        if event.mm_kin_vtx_prob[mm_index] < 0.1: continue
        if event.mm_kin_sl3d[mm_index] < 10: continue
        if event.mm_kin_slxy[mm_index] < 50: continue
        if event.mm_kin_pt[mm_index] > 15: continue
        # if event.mm_kin_alpha[mm_index]/event.mm_kin_alphaErr[mm_index] > 2: continue
        
        if event.mm_kin_lxy[mm_index] < 5: continue
        if event.MuonId_nLostHitsInner[event.mm_mu1_index[mm_index]] == 0: continue
        if event.MuonId_nLostHitsInner[event.mm_mu2_index[mm_index]] == 0: continue
        
        n_tight += 1

        # if event.mm_kin_lxy[mm_index] < 8: continue
        # if event.MuonId_nLostHitsInner[event.mm_mu1_index[mm_index]] < 2: continue
        # if event.MuonId_nLostHitsInner[event.mm_mu2_index[mm_index]] < 2: continue

        # n_very_tight += 1

        # if event.mm_kin_lxy[mm_index] < 20: continue
        
        # n_extreme += 1
        
        n_mus = 0
        if abs(event.mm_gen_mu1_pdgId[mm_index]) == 13: n_mus += 1
        if abs(event.mm_gen_mu2_pdgId[mm_index]) == 13: n_mus += 1

        if n_mus not in n_fakes:
            n_fakes[n_mus] = 0
        n_fakes[n_mus] += 1

        if verbose and (verbose_limit == None or verbose_count < verbose_limit):
            verbose_count += 1
            
            ancestor = event.mm_gen_cpdgId[mm_index]
            ancestor_index = event.mm_gen_cindex[mm_index]
            process = classify_event(event, ancestor, ancestor_index)
            if process not in processes:
                processes[process] = 0
            processes[process] += 1

            if print_only_unclassified and process != "Unclassified":
                continue
            
            print("mass: %0.3f, mass error: %0.3f, pt: %0.1f, alpha_signifance: %0.1f" % (
                event.mm_kin_mass[mm_index],
                event.mm_kin_massErr[mm_index],
                event.mm_kin_pt[mm_index],
                event.mm_kin_alpha[mm_index]/event.mm_kin_alphaErr[mm_index],
                )
            )
            
            if ancestor == 0:
                print("\tNo common ancestor")
            else:
                print("\tCommon ancestor: %d (event: %u, index: %d)" % (ancestor, event.event, ancestor_index))
                print("\t" + process)
                if ancestor_index >= 0:
                    id = ancestor_index
                    if origin_analysis:
                        # go up in ancestry
                        if event.GenPart_genPartIdxMother[id] >= 0:
                            id = event.GenPart_genPartIdxMother[id]
                        if event.GenPart_genPartIdxMother[id] >= 0:
                            id = event.GenPart_genPartIdxMother[id]
                    
                    for gen_index in range(event.nGenPart):
                        if event.GenPart_genPartIdxMother[gen_index] == id:
                            print_decay(event, gen_index, "\t\t")
                    
        # ## Determine origin
        # common_parrent = abs(event.mm_gen_cpdgId[mm_index])
        # if common_parrent not in signatures:
        #     signatures[common_parrent] = dict()
        
        # m1 = event.mm_gen_mu1_mpdgId[mm_index]
        # m2 = event.mm_gen_mu2_mpdgId[mm_index]

        # if abs(m1) > abs(m2):
        #     signature = name(m1) + " " + name(m2)
        # else:
        #     signature = name(m2) + " " + name(m1)

        # if signature not in signatures[common_parrent]:
        #     signatures[common_parrent][signature] = 0
        # signatures[common_parrent][signature] += 1

print("n_ref: ", n_ref)
print("n_tight: ", n_tight)
print("n_very_tight: ", n_very_tight)
print("n_extreme: ", n_extreme)

print("Event classification")
for process, n in sorted(processes.items(), key=lambda item: item[1], reverse=True):
    print(f"\t{process}: {n}")

# for common_parrent in sorted(signatures.keys(), key=lambda x:sum(signatures[x].values()), reverse=True): 
#     print("Common parrent: %s (N = %u)" % (common_parrent, sum(signatures[common_parrent].values())))
#     for signature in sorted(signatures[common_parrent], key=signatures[common_parrent].get, reverse=True):
#         print("\t%s: %u" % (signature, signatures[common_parrent][signature]))

# print("Unknown PDG ids:")
# print(unknown_pdgIds)

# print("Distribut by number of matched muons:")
# print(n_fakes)
