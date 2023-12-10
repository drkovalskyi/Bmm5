import ROOT
import glob
import os

# input_path = '/data/dmytro/Run3-Bmm-NanoAODv12/root-files/InclusiveDileptonMinBias.root'
input_path = '/tmp/dmytro/InclusiveDileptonMinBias.root'

verbose = True
force_recreate = False

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
    print("%spdgId: %d" % (indent, event.GenPart_pdgId[gen_index]))
    if event.GenPart_status[gen_index] != 1:
        for i in range(event.nGenPart):
            if event.GenPart_genPartIdxMother[i] == gen_index:
                print_decay(event, i, indent + "\t")


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

    df2 = df.Define("goodCandidates", "abs(mm_kin_mass - 0.5) < 0.02 && mm_kin_sl3d > 3 && mm_kin_alpha < 0.1")
    df_final = df2.Filter("Sum(goodCandidates) > 0", "Event has good candidates")

    df_final.Snapshot("Events", skim_file_name, "nMuon|Muon_.*|nmm|mm_.*|nGenPart|GenPart_.*|event|run")

signatures = dict()

n_fakes = dict()

f = ROOT.TFile(skim_file_name)
events = f.Get("Events")

if not events:
    raise Exception("Cannot read the skim")

print("Number of events in the skim: %d" % events.GetEntries())


# for event in chain:
for event in events:
    for mm_index in range(event.nmm):
        ## Selection
        if abs(event.mm_kin_mass[mm_index] - 0.5) > 0.02: continue
        if event.Muon_charge[event.mm_mu1_index[mm_index]] * \
           event.Muon_charge[event.mm_mu2_index[mm_index]] > 0: continue
        if event.Muon_softMva[event.mm_mu1_index[mm_index]] < 0.45: continue
        if event.Muon_softMva[event.mm_mu2_index[mm_index]] < 0.45: continue
        if event.mm_mu1_pt[mm_index] < 4: continue
        if event.mm_mu2_pt[mm_index] < 3: continue
        if event.mm_kin_vtx_prob[mm_index] < 0.01: continue
        if event.mm_kin_alpha[mm_index] > 0.01: continue
        if event.mm_kin_sl3d[mm_index] < 3: continue
        if event.mm_kin_lxy[mm_index] < 1: continue

        n_mus = 0
        if abs(event.mm_gen_mu1_pdgId[mm_index]) == 13: n_mus += 1
        if abs(event.mm_gen_mu2_pdgId[mm_index]) == 13: n_mus += 1

        if n_mus not in n_fakes:
            n_fakes[n_mus] = 0
        n_fakes[n_mus] += 1

        print("mass: %0.2f" % event.mm_kin_mass[mm_index])
        if verbose:
            if event.mm_gen_cpdgId[mm_index] == 0:
                print("\tNo common ancestor")
            else:
                ancestor_index = event.mm_gen_cindex[mm_index]
                print("\tCommon ancestor: %d (event: %u, index: %d)" % (event.mm_gen_cpdgId[mm_index], event.event, ancestor_index))
                if ancestor_index >= 0:
                    for gen_index in range(event.nGenPart):
                        if event.GenPart_genPartIdxMother[gen_index] == ancestor_index:
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

# for common_parrent in sorted(signatures.keys(), key=lambda x:sum(signatures[x].values()), reverse=True): 
#     print("Common parrent: %s (N = %u)" % (common_parrent, sum(signatures[common_parrent].values())))
#     for signature in sorted(signatures[common_parrent], key=signatures[common_parrent].get, reverse=True):
#         print("\t%s: %u" % (signature, signatures[common_parrent][signature]))

# print("Unknown PDG ids:")
# print(unknown_pdgIds)

# print("Distribut by number of matched muons:")
# print(n_fakes)
