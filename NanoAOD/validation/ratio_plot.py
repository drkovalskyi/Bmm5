import os
import ROOT
import math 
import sweights
import tdrstyle
from ROOT.RooFit import Binning
from ROOT import RooRealVar

# Set the TDR style
tdrstyle.setTDRStyle()

output_path = "/afs/cern.ch/user/d/dmytro/www/public_html/plots/AN/test/";


def print_canvas(output_name_without_extention, path, canvas=ROOT.gPad):
    if not os.path.exists(path):
        os.makedirs(path)
    canvas.Print("%s/%s.png" % (path, output_name_without_extention))
    canvas.Print("%s/%s.pdf" % (path, output_name_without_extention))
    canvas.Print("%s/%s.root" % (path, output_name_without_extention))
    canvas.Print("%s/%s.C" % (path, output_name_without_extention))


ROOT.gROOT.SetBatch(True)
c1 = ROOT.TCanvas("c2","c2", 800, 800)
# c2_1 = ROOT.TPad("c2_1", "", 0.0, 0.25, 1.0, 1.0)
# c2_2 = ROOT.TPad("c2_2", "", 0.0, 0.0, 1.0, 0.25)

chain = ROOT.TChain("mva")
chain.Add("/eos/cms/store/group/phys_bphys/bmm/bmm5/PostProcessing/FlatNtuples/516/bmm_mva_jpsik/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v2+MINIAODSIM/0beb1bfffb548eeef0bbd03a4eecf864.root")

h1 = ROOT.TH1F("h1", ";p_{T}(#mu#mu)", 30, 0, 30)
h2 = ROOT.TH1F("h2", "", 30, 0, 30)
chain.Draw("mm_kin_pt>>h1", "mm_mva>0.9")
chain.Draw("mm_kin_pt>>h2", "mm_mva<0.9")

h1.SetLineWidth(2)
h1.SetMarkerStyle(21)

            
h2.SetMarkerColor(ROOT.kBlue)
h2.SetMarkerStyle(20)

h2.Scale(h1.Integral()/h2.Integral())
max_value = h1.GetMaximum()
if h2.GetMaximum() > max_value:
    max_value = h2.GetMaximum()
h1.SetMaximum(max_value * 1.1)
h1.SetMinimum(0)
h2.SetMaximum(max_value * 1.1)
# c2_1.cd()
h1.Draw("hist e")
h2.Draw("same")

legend = ROOT.TLegend(0.15,0.75,0.5,0.87)
legend.SetFillStyle(0)
legend.SetLineWidth(0)
legend.SetBorderSize(1)

legend.AddEntry(h1, "MVA>0.9")
legend.AddEntry(h2, "MVA<0.9")
legend.Draw()


ratio_plot = ROOT.TRatioPlot(h1, h2)
ratio_plot.SetH1DrawOpt("hist e");
ratio_plot.Draw()
# ratio_plot.GetUpperPad().BuildLegend(0.15,0.75,0.5,0.87)

ratio_plot.SetSeparationMargin(0.03)
ratio_plot.GetLowerRefGraph().SetMinimum(0.4)
ratio_plot.GetLowerRefGraph().SetMaximum(1.6)
# ratio_plot.GetXaxis().SetTitleSize()
# SetBottomMargin(2.0)
# c1.SetBottomMargin(2.0)
# rp->GetLowerRefYaxis()->SetRange(...)
ratio_plot.GetLowerRefYaxis().SetTitle("Data/MC")

ratio_plot.GetLowerRefYaxis().SetTitleSize()
ratio_plot.GetLowerRefYaxis().SetTitleOffset(1.1)
ratio_plot.GetLowerRefYaxis().SetLabelSize(0.035)
ratio_plot.GetLowYaxis().SetNdivisions(503)

ratio_plot.GetLowerRefXaxis().SetTitleSize()
ratio_plot.GetLowerRefXaxis().SetTitleOffset()
ratio_plot.GetLowerRefXaxis().SetLabelSize(0.035)

ratio_plot.GetUpperRefYaxis().SetTitleSize()
ratio_plot.GetUpperRefYaxis().SetTitleOffset()
ratio_plot.GetUpperRefYaxis().SetLabelSize(0.035)

ratio_plot.GetUpperRefXaxis().SetTitleSize()
ratio_plot.GetUpperRefXaxis().SetTitleOffset(0)
ratio_plot.GetUpperRefXaxis().SetLabelSize(0.035)

c1.cd()
legend.Draw()

# rp1.GetUpperRefYaxis().SetTitle("entries")
c1.Update()

#     # print_canvas(var.GetName() + "_ratio_" + name, output_path)


# c2_2.cd()
# h_ratio = h1.Clone("h_ratio")
# h_ratio.Divide(h1, h2)
# h_ratio.SetMinimum(0.)
# h_ratio.SetMaximum(2.0)
# h_ratio.Draw()
# c2.cd()
# c2_1.Draw()
# c2_2.Draw()

print_canvas("test", output_path)

#     # print_canvas(name + "_ratio_" + v.GetName(), output_path)
#     print_canvas(var.GetName() + "_ratio_" + name, output_path)

