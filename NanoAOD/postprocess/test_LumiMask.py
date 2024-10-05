import ROOT
import json

# Function to convert JSON to custom string format
def json_to_custom_format(filename):
    try:
        with open(filename, 'r') as file:
            json_data = json.load(file)
        
        parts = []
        for run, lumi_ranges in json_data.items():
            lumi_parts = [f"{start}-{end}" for start, end in lumi_ranges]
            parts.append(f"{run}:{','.join(lumi_parts)}")
        return ';'.join(parts)
    except json.JSONDecodeError as e:
        print(f"JSON parsing error: {e}")
        return None

# Load lumi masks
lumi_mask_string = ""
for json_filename in [
        "certification/muon/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON_MuonPhys.txt",
        "certification/muon/Cert_294927-306462_13TeV_UL2017_Collisions17_MuonJSON.txt",
        "certification/muon/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON_MuonPhys.txt",
        "certification/muon/Cert_Collisions2022_355100_362760_Muon.json",
        "certification/muon/Cert_Collisions2023_366442_370790_Muon.json",
]:
    if lumi_mask_string != "":
        lumi_mask_string += ';'
    lumi_mask_string += json_to_custom_format(json_filename)
# print(lumi_mask_string)

# Check if conversion succeeded
if lumi_mask_string:
    # Declare the LumiMask class in ROOT using the header file
    with open('LumiMask.h', 'r') as file:
        lumi_mask_code = file.read()
    ROOT.gInterpreter.Declare(lumi_mask_code)

    # Create the LumiMask object using the custom string format
    lumi_mask = ROOT.LumiMask.fromCustomString(lumi_mask_string, 0, 0)

    # Perform some test cases
    test_cases = [
        (355100, 100),  # Expected: True if 100 is within the range in JSON
        (355101, 5),    # Expected: True if 5 is within the range in JSON
        (355102, 30),   # Expected: False if 30 is outside the range in JSON
        (355103, 1)     # Expected: True if 1 is within the range in JSON
    ]

    for run, lumi in test_cases:
        if lumi_mask.accept(run, lumi):
            print(f"LumiBlock {lumi} in Run {run} is accepted.")
        else:
            print(f"LumiBlock {lumi} in Run {run} is not accepted.")
else:
    print("Failed to convert JSON to custom string format.")

# Test RDataFrame
tree = ROOT.TChain("Events")
tree.Add("/eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/529/ParkingDoubleMuonLowMass0+Run2022C-PromptReco-v1+MINIAOD/7f06673d-d50c-4bdc-b37c-d75d0ca5565d.root")
tree.Add("/eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/529/HLTPhysics+Run2022C-PromptReco-v1+MINIAOD/8107f6dc-4bb0-4039-b864-9c014fa6aaa9.root")

ROOT.gInterpreter.Declare(f'''
std::string lumi_mask_string;

bool passed_lumi_mask(unsigned int run, unsigned int lumi) {{
    static LumiMask lumi_mask = LumiMask::fromCustomString(lumi_mask_string);
    return lumi_mask.accept(run, lumi);
}}
''')

ROOT.gInterpreter.ProcessLine(f'lumi_mask_string = "{lumi_mask_string}";')


df = ROOT.RDataFrame(tree)
n_events = df.Count().GetValue()
print("Number of events to process: %d" % n_events)

df = df.Define("certified", "passed_lumi_mask(run, luminosityBlock)")
df = df.Filter("certified == 1", "passed data certification")
n_passed = df.Count().GetValue()
print(f"Number of events passing the golden JSON: {n_passed}")
