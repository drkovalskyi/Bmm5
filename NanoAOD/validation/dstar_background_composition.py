import ROOT
import glob

input_path = '/eos/cms/store/group/phys_bphys/bmm/bmm6/PostProcessing/Skims/526/dstar_mm/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8+Run3Summer22MiniAODv3-Pilot_124X_mcRun3_2022_realistic_v12-v5+MINIAODSIM/'

pdgIds = {
    511:'B0', 521:'B+', 531:'Bs', 541:'Bc',
    411:'D0', 421:'D+', 431:'Ds',
    321:'K+', 211:'pi+', 15:'tau',
}

unknown_pdgIds = []

def name(pdgId):
    id = abs(pdgId)
    if id not in pdgIds:
        if id not in unknown_pdgIds:
            unknown_pdgIds.append(id)
        return "other"
    return pdgIds[id]


chain = ROOT.TChain("Events")
for f in glob.glob("%s/*.root" % input_path):
    chain.Add(f)

signatures = dict()

n_fakes = dict()

for event in chain:
    for dstar_index in range(event.ndstar):
        ## Selection
        if event.dstar_mm_index[dstar_index] < 0: continue
        mm_index = event.dstar_mm_index[dstar_index]
        if event.mm_kin_alpha[mm_index] > 0.1: continue
        if event.mm_kin_sl3d[mm_index] < 3: continue
        if event.mm_kin_vtx_prob[mm_index] < 0.01: continue
        if event.dstar_pv_with_pion_prob[dstar_index] < 0.1: continue
        if event.dstar_dm_pv[dstar_index] < 0.140 or event.dstar_dm_pv[dstar_index] > 0.155: continue
        if event.mm_kin_mass[mm_index] < 1.76 or event.mm_kin_mass[mm_index] > 1.94: continue

        n_mus = 0
        if abs(event.mm_gen_mu1_pdgId[mm_index]) == 13: n_mus += 1
        if abs(event.mm_gen_mu2_pdgId[mm_index]) == 13: n_mus += 1

        if n_mus not in n_fakes:
            n_fakes[n_mus] = 0
        n_fakes[n_mus] += 1

        ## Classify only events were both muons are matched
        if n_mus < 2: continue
        
        ## Determine origin
        common_parrent = abs(event.mm_gen_cpdgId[mm_index])
        if common_parrent not in signatures:
            signatures[common_parrent] = dict()
        
        m1 = event.mm_gen_mu1_mpdgId[mm_index]
        m2 = event.mm_gen_mu2_mpdgId[mm_index]

        if abs(m1) > abs(m2):
            signature = name(m1) + " " + name(m2)
        else:
            signature = name(m2) + " " + name(m1)

        if signature not in signatures[common_parrent]:
            signatures[common_parrent][signature] = 0
        signatures[common_parrent][signature] += 1

for common_parrent in sorted(signatures.keys(), key=lambda x:sum(signatures[x].values()), reverse=True): 
    print("Common parrent: %s (N = %u)" % (common_parrent, sum(signatures[common_parrent].values())))
    for signature in sorted(signatures[common_parrent], key=signatures[common_parrent].get, reverse=True):
        print("\t%s: %u" % (signature, signatures[common_parrent][signature]))

print("Unknown PDG ids:")
print(unknown_pdgIds)

print("Distribut by number of matched muons:")
print(n_fakes)
