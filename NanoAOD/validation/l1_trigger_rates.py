import ROOT
import re

# files = [
#     '/eos/cms/store/data/Run2023D/ZeroBias/NANOAOD/PromptReco-v1/000/370/580/00000/127df71f-8c3b-4d1b-88f0-b8cc983c5d3b.root',
#     '/eos/cms/store/data/Run2023D/ZeroBias/NANOAOD/PromptReco-v1/000/370/580/00000/8489c431-473e-4378-a2fe-9c8afe703562.root',
#     '/eos/cms/store/data/Run2023D/ZeroBias/NANOAOD/PromptReco-v1/000/370/580/00000/04fdb3c3-5460-4b6e-826b-76407d5674db.root',
#     '/eos/cms/store/data/Run2023D/ZeroBias/NANOAOD/PromptReco-v1/000/370/580/00000/ef893e7b-ece2-4bf5-86bd-cbf66a32d80f.root',
#     '/eos/cms/store/data/Run2023D/ZeroBias/NANOAOD/PromptReco-v1/000/370/580/00000/17e6c1b0-9aaf-4d1b-83ab-dbc7cf006138.root',
#     '/eos/cms/store/data/Run2023D/ZeroBias/NANOAOD/PromptReco-v1/000/370/580/00000/637a67a1-54b8-4dfd-a92f-a98af78a37b3.root',
#     '/eos/cms/store/data/Run2023D/ZeroBias/NANOAOD/PromptReco-v1/000/370/580/00000/80fbb1d8-22e1-4485-b597-f532ca6b2605.root',
#     '/eos/cms/store/data/Run2023D/ZeroBias/NANOAOD/PromptReco-v1/000/370/580/00000/50e1442c-a7f5-4cd4-97bf-ad05dd95afb9.root',
#     '/eos/cms/store/data/Run2023D/ZeroBias/NANOAOD/PromptReco-v1/000/370/580/00000/c5de979c-7ec3-46ac-b3bf-7d4c725bb362.root',
#     '/eos/cms/store/data/Run2023D/ZeroBias/NANOAOD/PromptReco-v1/000/370/580/00000/14773b1f-f167-4489-aae3-d428f0c71e46.root',
# ]
# scale_factor = 29.5 # Hz / event_count
# filter = "HLT_ZeroBias"

files = [
    '/eos/cms/store/data/Run2023D/ZeroBias/NANOAOD/PromptReco-v2/000/370/772/00000/c34d019e-e702-4e7f-acfe-421cf9aac8a4.root',
    '/eos/cms/store/data/Run2023D/ZeroBias/NANOAOD/PromptReco-v2/000/370/772/00000/7e2dd65f-ba87-40f7-a84a-f878377a8978.root',
    '/eos/cms/store/data/Run2023D/ZeroBias/NANOAOD/PromptReco-v2/000/370/772/00000/88c3bd03-f25b-420c-96a5-b9f554212ee0.root',
    '/eos/cms/store/data/Run2023D/ZeroBias/NANOAOD/PromptReco-v2/000/370/772/00000/b22ea534-d353-4314-9978-80b496a39e9b.root',
    '/eos/cms/store/data/Run2023D/ZeroBias/NANOAOD/PromptReco-v2/000/370/772/00000/0cdb595d-139d-4a8c-b977-b10d3690e474.root',
]
scale_factor = 49.9 # Hz / event_count
filter = "HLT_ZeroBias && luminosityBlock >= 29 && luminosityBlock <= 1247"

masked_l1s = [
    'L1_ZeroBias',
    'L1_UnprefireableEvent',
    # double electron
    'L1_DoubleEG4_er1p2_dR_Max0p9',
    'L1_DoubleEG4p5_er1p2_dR_Max0p9',
    'L1_DoubleEG5_er1p2_dR_Max0p9',
    'L1_DoubleEG5p5_er1p2_dR_Max0p8',
    'L1_DoubleEG6_er1p2_dR_Max0p8',
    'L1_DoubleEG6p5_er1p2_dR_Max0p8',
    'L1_DoubleEG7_er1p2_dR_Max0p8',
    'L1_DoubleEG7p5_er1p2_dR_Max0p7',
    'L1_DoubleEG8_er1p2_dR_Max0p7',
    'L1_DoubleEG10_er1p2_dR_Max0p6',
    'L1_DoubleEG10p5_er1p2_dR_Max0p6',
    'L1_DoubleEG8p5_er1p2_dR_Max0p7',
    'L1_DoubleEG9_er1p2_dR_Max0p7',
    'L1_DoubleEG9p5_er1p2_dR_Max0p6',
]

bph_2023_l1s = [
    ## double muon
    'L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4',
    'L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6',
    'L1_DoubleMu4p5_SQ_OS_dR_Max1p2',
    'L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4',
    'L1_DoubleMu3er2p0_SQ_OS_dR_Max1p6',
    'L1_DoubleMu4_SQ_OS_dR_Max1p2',
    ## tripple muon
    'L1_TripleMu_3SQ_2p5SQ_0_OS_Mass_Max12',
    'L1_TripleMu_4SQ_2p5SQ_0_OS_Mass_Max12',
    'L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9',
    'L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9',
    ## double electron
    # prescaled 
    # 'L1_DoubleEG4_er1p2_dR_Max0p9',
    # 'L1_DoubleEG4p5_er1p2_dR_Max0p9',
    # 'L1_DoubleEG5_er1p2_dR_Max0p9',
    # 'L1_DoubleEG5p5_er1p2_dR_Max0p8',
    # 'L1_DoubleEG6_er1p2_dR_Max0p8',
    # 'L1_DoubleEG6p5_er1p2_dR_Max0p8',
    # 'L1_DoubleEG7_er1p2_dR_Max0p8',
    # 'L1_DoubleEG7p5_er1p2_dR_Max0p7',
    # 'L1_DoubleEG8_er1p2_dR_Max0p7',
    # 'L1_DoubleEG10_er1p2_dR_Max0p6',
    # 'L1_DoubleEG10p5_er1p2_dR_Max0p6',
    # 'L1_DoubleEG8p5_er1p2_dR_Max0p7',
    # 'L1_DoubleEG9_er1p2_dR_Max0p7',
    # 'L1_DoubleEG9p5_er1p2_dR_Max0p6',
    # unprescaled 
    'L1_DoubleEG11_er1p2_dR_Max0p6',
]

l1s_to_drop = [
    ## double muon
    'L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4',
    'L1_DoubleMu3er2p0_SQ_OS_dR_Max1p6',
    'L1_DoubleMu4_SQ_OS_dR_Max1p2',
    ## tripple muon
    'L1_TripleMu_3SQ_2p5SQ_0_OS_Mass_Max12',
    'L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9',
    ## double electron
    # prescaled 
    # 'L1_DoubleEG4_er1p2_dR_Max0p9',
    # 'L1_DoubleEG4p5_er1p2_dR_Max0p9',
    # 'L1_DoubleEG5_er1p2_dR_Max0p9',
    # 'L1_DoubleEG5p5_er1p2_dR_Max0p8',
    # 'L1_DoubleEG6_er1p2_dR_Max0p8',
    # 'L1_DoubleEG6p5_er1p2_dR_Max0p8',
    # 'L1_DoubleEG7_er1p2_dR_Max0p8',
    # 'L1_DoubleEG7p5_er1p2_dR_Max0p7',
    # 'L1_DoubleEG8_er1p2_dR_Max0p7',
    # 'L1_DoubleEG10_er1p2_dR_Max0p6',
    # 'L1_DoubleEG10p5_er1p2_dR_Max0p6',
    # 'L1_DoubleEG8p5_er1p2_dR_Max0p7',
    # 'L1_DoubleEG9_er1p2_dR_Max0p7',
    # 'L1_DoubleEG9p5_er1p2_dR_Max0p6',
    # unprescaled 
    'L1_DoubleEG11_er1p2_dR_Max0p6',
    ## single EG
    # HLT_Ele30_WPTight_Gsf_v5 -> HLT_Ele32_WPTight_Gsf_v19
    'L1_SingleIsoEG30er2p5',
    'L1_SingleIsoEG30er2p1',
    # other
    'L1_SingleEG36er2p5'
]


all_l1s = []

chain = ROOT.TChain("Events")
for f in files:
    chain.Add(f)

for br in chain.GetListOfBranches():
    name = br.GetName()
    if re.search(r"^L1_", name):
        # print(name)
        if name not in masked_l1s: 
            all_l1s.append(name)

# print(all_l1s)

df = ROOT.RDataFrame(chain)
df_zb = df.Filter(filter)
n_events = df_zb.Count().GetValue()
    
# print("Number of certified events passed HLT_ZeroBias trigger: %u, \trate: %0.1f kHz" % (n_events, n_events * scale_factor / 1e3))

bph_2024_l1s = []

for l1 in sorted(bph_2023_l1s):
    print("%s rate: %0.2f kHz" % (l1, df_zb.Filter(l1).Count().GetValue() * scale_factor / 1e3))
    if l1 not in l1s_to_drop:
        bph_2024_l1s.append(l1)

bph_2023_l1_rate = df_zb.Filter(" || ".join(bph_2023_l1s)).Count().GetValue() * scale_factor / 1e3
print("BPH 2023 L1 rate: %0.2f kHz" % (bph_2023_l1_rate))

bph_2024_l1_rate = df_zb.Filter(" || ".join(bph_2024_l1s)).Count().GetValue() * scale_factor / 1e3
print("BPH 2024 L1 rate: %0.2f kHz" % (bph_2024_l1_rate))

print("Saved rate: %0.1f kHz" % (bph_2023_l1_rate - bph_2024_l1_rate))

###############################
## Full menu
###############################

# for l1 in sorted(all_l1s):
#     print("%s rate: %0.1f kHz" % (l1, df_zb.Filter(l1).Count().GetValue() * scale_factor / 1e3))

l1_Run2023_rate = df_zb.Filter(" || ".join(all_l1s)).Count().GetValue() * scale_factor / 1e3
print("L1 menu rate: %0.1f kHz" % l1_Run2023_rate)

l1_new_menu = []
for l1 in all_l1s:
    if l1 in l1s_to_drop:
        continue
    l1_new_menu.append(l1)

l1_new_menu_rate = df_zb.Filter(" || ".join(l1_new_menu)).Count().GetValue() * scale_factor / 1e3
print("New L1 menu rate: %0.1f kHz" % l1_new_menu_rate)
print("Saved rate: %0.1f kHz" % (l1_Run2023_rate - l1_new_menu_rate))
