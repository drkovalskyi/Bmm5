import ROOT, glob, os, json

prefix = "/eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/526/"

input_paths = [
    # prefix + "ZeroBias+Run2023C-PromptReco-v1+MINIAOD",
    # prefix + "ZeroBias+Run2023C-PromptReco-v2+MINIAOD",
    # prefix + "ZeroBias+Run2023C-PromptReco-v3+MINIAOD",
    # prefix + "ZeroBias+Run2023C-PromptReco-v4+MINIAOD",
    prefix + "ZeroBias+Run2023D-PromptReco-v1+MINIAOD",
    prefix + "ZeroBias+Run2023D-PromptReco-v2+MINIAOD",
]

trigger = "HLT_ZeroBias" 

goodruns = dict()

def is_certified(event, type):
    # load certification information
    if type not in goodruns:
        goodruns[type] = dict()
        for f in os.listdir('certification/%s' % type):
            goodruns[type].update(json.load(open('certification/%s/%s' % (type, f))))
        print("Number of runs in the %s certification: %u" % (type, len(goodruns[type])))

    # run number is a string for some reason
    run = str(event.run)

    if run not in goodruns[type]:
        return False

    for min_lumi, max_lumi in goodruns[type][run]:
        if event.luminosityBlock >= min_lumi and event.luminosityBlock <= max_lumi:
            return True
    return False

chain = ROOT.TChain("Events")
for path in input_paths:
    for f in glob.glob("%s/*.root" % path):
        chain.Add(f)

print("n_total: %u" % chain.GetEntries())

# make a skim
df = ROOT.RDataFrame(chain)
df_trigger = df.Filter(trigger)

df_trigger.Snapshot("Events", "/tmp/dmytro/tmp_tree.root", ["run", "luminosityBlock"]);

f = ROOT.TFile("/tmp/dmytro/tmp_tree.root")
tree = f.Get("Events")
print("n_triggered: %u" % tree.GetEntries())

n_certified = 0
for event in tree:
    if not is_certified(event, "muon"):
        continue
    n_certified += 1

print("n_certified: %u" % n_certified)


