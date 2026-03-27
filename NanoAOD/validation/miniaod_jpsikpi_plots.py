import sys,os,re
import ROOT
from DataFormats.FWLite import Events, Handle
from math import *
import numpy

events = Events (['/eos/cms/store/user/dmytro/tmp/store+mc+RunIISummer20UL18MiniAODv2+BdToJPsiKPi_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+MINIAODSIM+106X_upgrade2018_realistic_v16_L1v1-v3+2820000+4672866D-5E71-2E46-9A0E-4CB8EE1ED3D9.root'])
nEvents = events.size()
print("Number of events: %d" % nEvents)

output_path = "/eos/home-d/dmytro/www/plots/tmp/jspikpi/"

def print_canvas(output_name_without_extention, path, canvas=ROOT.gPad):
    if not os.path.exists(path):
        os.makedirs(path)
    canvas.Print(f"{path}/{output_name_without_extention}.png")
    canvas.Print(f"{path}/{output_name_without_extention}.pdf")
    canvas.Print(f"{path}/{output_name_without_extention}.root")


handlePruned  = Handle ("std::vector<reco::GenParticle>")
labelPruned = ("prunedGenParticles")

h_kpi_mass_all = ROOT.TH1D("h_kpi_mass_all","Gen Kpi mass; mass, GeV", 100, 0, 2.5)
h_kpi_mass_1   = ROOT.TH1D("h_kpi_mass_1","Gen Kpi mass (track pt>1GeV);mass, GeV", 100, 0, 2.5)
h_kpi_mass_1z  = ROOT.TH1D("h_kpi_mass_1z","Gen Kpi mass (track pt>1GeV);mass, GeV", 25, 0.816, 0.976)
i = 0
for event in events:
    i += 1
    if i%10000==0:
        sys.stdout.write(".")
        sys.stdout.flush()
    # if i>30: break
    event.getByLabel (labelPruned, handlePruned)
    pruned = handlePruned.product()
    jpsi = None
    kaon = None
    pion = None

    # https://github.com/cms-sw/cmssw/blob/CMSSW_7_4_X/DataFormats/HepMCCandidate/interface/GenStatusFlags.h
    for p in pruned :
        if abs(p.pdgId()) == 443 and abs(p.mother().pdgId()) == 511:
            jpsi = p
        if abs(p.pdgId()) == 321 and abs(p.mother().pdgId()) == 511:
            kaon = p
        if abs(p.pdgId()) == 211 and abs(p.mother().pdgId()) == 511:
            pion = p
            
    # if jpsi and kaon and pion and kaon.charge() != pion.charge():
    if jpsi.mother() == kaon.mother() == pion.mother():
        kpi_mass = (kaon.p4() + pion.p4()).mass()
        h_kpi_mass_all.Fill(kpi_mass)
        if kaon.pt() > 1 and pion.pt() > 1:
            h_kpi_mass_1.Fill(kpi_mass)
            h_kpi_mass_1z.Fill(kpi_mass)
        # print(f"kpi_mass: {kpi_mass:0.1f}")
sys.stdout.write("\n")
sys.stdout.flush()

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetLabelSize(0.045, "X")
ROOT.gStyle.SetLabelSize(0.045, "Y")
ROOT.gStyle.SetTitleSize(0.045, "X")
ROOT.gStyle.SetTitleSize(0.045, "Y")
ROOT.gStyle.SetTitleOffset(1.2, "X")
ROOT.gStyle.SetTitleOffset(1.2, "Y")
ROOT.gStyle.SetPadLeftMargin(0.15)
ROOT.gStyle.SetPadRightMargin(0.15)
ROOT.gStyle.SetPadBottomMargin(0.15)
ROOT.gStyle.SetOptFit(1)

c = ROOT.TCanvas("c", "", 600, 600)

h_kpi_mass_all.Draw()
print_canvas("h_kpi_mass_all", output_path)
h_kpi_mass_1.Draw()
print_canvas("h_kpi_mass_1", output_path)
print(f"N events (h_kpi_mass_1): {h_kpi_mass_1.Integral():0.1f}")
h_kpi_mass_1z.Draw()
h_kpi_mass_1z.SetMinimum(0)
print(f"N events (h_kpi_mass_1z): {h_kpi_mass_1z.Integral():0.1f}")
print_canvas("h_kpi_mass_1z", output_path)
