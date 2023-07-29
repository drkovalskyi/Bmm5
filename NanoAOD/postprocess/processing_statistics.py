import ROOT, copy, subprocess, re

path = "/eos/cms/store/group/phys_bphys/bmm/bmm5/PostProcessing/FlatNtuples/518"

datasets = {}

command = 'find -L %s -maxdepth 2 -type d' % (path)
all_folders = subprocess.check_output(command, shell=True, encoding='utf8').splitlines()

pattern = 'Charmonium'

for folder in all_folders:
    if pattern != "":
        if not re.search(pattern, folder):
            continue
    match = re.search("([^\/]+)\/([^\/]+)\/*$", folder)
    if match:
        type = match.group(1)
        ds = match.group(2)

        if not ds in datasets:
            datasets[ds] = []
        datasets[ds].append(type) 

for ds in sorted(datasets):
    print(ds)
    for type in sorted(datasets[ds]):
        print("\t%-30s " % type, end=' ')
        chain = ROOT.TChain("info")
        chain.Add("%s/%s/%s/*.root" % (path, type, ds))
        n_processed = 0
        for entry in chain:
            n_processed += entry.n_processed
        print(" \t:", n_processed)
    
