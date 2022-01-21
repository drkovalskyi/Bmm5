import os
import ROOT
import tdrstyle

# Set the TDR style
tdrstyle.setTDRStyle()

version = 516

output_path = "/afs/cern.ch/user/d/dmytro/www/public_html/plots/bmm5_NanoAODv8-%u/random/" % version

path = "/eos/cms/store/group/phys_bphys/bmm/bmm5/NanoAOD/%u/" % version



def print_canvas(output_name_without_extention, path, canvas=ROOT.gPad):
    if not os.path.exists(path):
        os.makedirs(path)
    canvas.Print("%s/%s.png" % (path, output_name_without_extention))
    canvas.Print("%s/%s.pdf" % (path, output_name_without_extention))
    canvas.Print("%s/%s.root" % (path, output_name_without_extention))
    # canvas.Print("%s/%s.C" % (path, output_name_without_extention))

ROOT.gROOT.SetBatch(True)
c1 = ROOT.TCanvas("c1","c1", 800, 600)

name = "bhh_mass"
selection = "mm_mu1_index>=0 && mm_mu2_index>=0 && abs(mm_kin_mu1eta)<1.4 && mm_kin_mu1pt>4" +\
    " && abs(mm_kin_mu2eta)<1.4 && mm_kin_mu2pt>4 && mm_kin_sl3d>4 && mm_kin_pt>5.0" +\
    " && mm_kin_vtx_prob>0.025"

chain = ROOT.TChain("Events")
print "Number of files:",
print chain.Add(path + "/BTohh_hToMuNu_BsBdMixture_modHadLifetime_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v3+MINIAODSIM/*.root")
h = ROOT.TH1F("h",";m_{#mu#mu}, [GeV]", 100, 4.9, 5.9)
chain.Draw("mm_kin_mass>>h", selection)
h.SetLineWidth(2)
h.SetFillColor(ROOT.kMagenta)
h.Draw()

print_canvas(name, output_path)
