import os
import ROOT
import tdrstyle

tdrstyle.setTDRStyle()
ROOT.gROOT.SetBatch(True)
c1 = ROOT.TCanvas("c1","c1", 1200, 800)

selection = "mm_mu1_index>=0 && mm_mu2_index>=0 && "\
    "Muon_softMva[mm_mu1_index]>0.45 && "\
    "abs(mm_kin_mu1eta)<1.4 && "\
    "mm_kin_mu1pt>4 && "\
    "Muon_softMva[mm_mu2_index]>0.45 && "\
    "abs(mm_kin_mu2eta)<1.4 && "\
    "mm_kin_mu2pt>4 && "\
    "abs(mm_kin_mass-%f)<0.5 && "\
    "(abs(mm_kin_mu1eta)>0.7||abs(mm_kin_mu2eta)>0.7)"
    # "mm_kin_sl3d>6 && "\
    # "mm_kin_pt>5.0 && mm_kin_vtx_prob>0.025 && "\

# selection += " && mm_mva>0.99"
# selection += " && mm_gen_pdgId!=0"

input_path = "/eos/cms/store/group/phys_bphys/bmm/bmm5/NanoAOD/517/"
output_path = "/afs/cern.ch/user/d/dmytro/www/public_html/plots/Bmm/AN/tmp/"

def print_canvas(output_name_without_extention, path, canvas=ROOT.gPad):
    if not os.path.exists(path):
        os.makedirs(path)
    canvas.Print("%s/%s.png" % (path, output_name_without_extention))
    canvas.Print("%s/%s.pdf" % (path, output_name_without_extention))
    # canvas.Print("%s/%s.root" % (path, output_name_without_extention))
    # canvas.Print("%s/%s.C" % (path, output_name_without_extention))

samples = [
    # {
    #     'name':'Bs',
    #     'files': input_path + '/BsToMuMu_SoftQCDnonD_TuneCP5_BsLifetime1p40_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v1+MINIAODSIM/0007E4A0-4C83-0140-8EC3-973B13CD9F27.root'
    # },
    {
        'name':'Jpsi',
        'mass':3.097,
        'files': input_path + '/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v2+MINIAODSIM/0*.root'
    },
]

def process(sample):
    print "Processing sample: %s" % sample['name']
    
    events = ROOT.TChain("Events")
    print "Number of files:", events.Add(sample['files'])

    prof = ROOT.TProfile("prof","", 50, 0, 50)
    prof.SetMarkerStyle(20)
    prof.GetYaxis().SetRangeUser(0.980,1.010)

    events.Draw("mm_kin_mass/%f:mm_kin_pt>>prof" % sample['mass'], selection % sample['mass'])
    print_canvas("mc_mass_scale_%s" % sample['name'], output_path)

for sample in samples:
    process(sample)


#     jpi_bf = 3.92e-5
# # jpi_f = ROOT.TFile.Open("/eos/cms/store/group/phys_bphys/bmm/bmm5/NanoAOD/516/BuToJpsiPi_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v1+MINIAODSIM/00956A81-A197-6A4A-AF29-701466920A11.root")
# # jpi_lumis = jpi_f.Get("LuminosityBlocks")
# # jpi_events = jpi_f.Get("Events")
# jpi_lumis = ROOT.TChain("LuminosityBlocks")
# jpi_lumis.Add("/eos/cms/store/group/phys_bphys/bmm/bmm5/NanoAOD/516/BuToJpsiPi_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v1+MINIAODSIM/*.root")
# jpi_events = ROOT.TChain("Events")
# jpi_events.Add("/eos/cms/store/group/phys_bphys/bmm/bmm5/NanoAOD/516/BuToJpsiPi_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v1+MINIAODSIM/*.root")

# n_jpi = 0
# for lumi in jpi_lumis:
#     n_jpi += lumi.GenFilter_numEventsTotal

# print n_jpi

# jpi_h = ROOT.TH1F("jpi_h","", 100, 4.9, 5.9)
# jpi_h.SetLineWidth(2)
# jpi_h.SetFillColor(ROOT.kYellow)
# jpi_events.Draw("bkmm_jpsimc_mass>>jpi_h", selection)

# jpi_h2 = ROOT.TH1F("jpi_h2","", 50, 0.95, 1.0)
# jpi_h2.SetLineWidth(2)
# jpi_h2.SetLineColor(ROOT.kBlue)
# jpi_h2.GetXaxis().SetNdivisions(505)
# jpi_events.Draw("bkmm_bmm_mva>>jpi_h2", selection)

# print_canvas("hist_bjpsipi_mva", output_path)

# jk_bf = 1.02e-3
# # jk_f = ROOT.TFile.Open("/eos/cms/store/group/phys_bphys/bmm/bmm5/NanoAOD/516/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v2+MINIAODSIM/01033CC8-1059-E54A-AE50-913E0E9BBF45.root")
# # jk_lumis = jk_f.Get("LuminosityBlocks")
# # jk_events = jk_f.Get("Events")
# jk_lumis = ROOT.TChain("LuminosityBlocks")
# jk_lumis.Add("/eos/cms/store/group/phys_bphys/bmm/bmm5/NanoAOD/516/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v2+MINIAODSIM/*.root")
# jk_events = ROOT.TChain("Events")
# jk_events.Add("/eos/cms/store/group/phys_bphys/bmm/bmm5/NanoAOD/516/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v2+MINIAODSIM/*.root")

# n_jk = 0
# for lumi in jk_lumis:
#     n_jk += lumi.GenFilter_numEventsTotal

# print n_jk

# jk_h = ROOT.TH1F("jk_h","", 100, 4.9, 5.9)
# jk_h.SetLineWidth(2)
# jk_h.SetFillColor(ROOT.kBlue)
# jk_events.Draw("bkmm_jpsimc_mass>>jk_h", selection)

# print_canvas("hist_bjpsik", output_path)

# jk_h2 = ROOT.TH1F("jk_h2","", 50, 0.95, 1.0)
# jk_h2.SetLineWidth(2)
# jk_h2.GetXaxis().SetNdivisions(505)
# jk_events.Draw("bkmm_bmm_mva>>jk_h2", selection)

# print_canvas("hist_bjpsik_mva", output_path)

# # stack
# scale = float(n_jk)/n_jpi*jpi_bf/jk_bf
# jpi_h.Scale(scale)

# hs = ROOT.THStack("hs","");
# hs.Add(jk_h)
# hs.Add(jpi_h)
# hs.Draw("hist")

# print_canvas("hist_bjpsiX", output_path)

# jk_h2.Draw()
# jpi_h2.Scale(jk_h2.Integral()/jpi_h2.Integral())
# jpi_h2.Draw("same e")

