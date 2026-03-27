import ROOT
import glob
from collections import defaultdict
import os

input_path = '/eos/cms/store/group/phys_bphys/bmm/bmm6/PostProcessing/FlatNtuples/532/fit-bkstarmm/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8+Run3Summer22MiniAODv4-validDigi_130X_mcRun3_2022_realistic_v5-v4+MINIAODSIM/'
output_path = "/eos/home-d/dmytro/www/plots/tmp/jspikpi/"


pdgIds = {
    511:'B0', 521:'B+', 531:'Bs', 541:'Bc',
    411:'D0', 421:'D+', 431:'Ds', 443:'Jpsi',
    10313:'K1(1270)', 10323:'K1(1270)+',
    321:'K+', 313:'K*', 333:'phi', 323:'K*+',
    315:'K*2(1430)', 331:"eta'", 221:'eta',
    2212:'p', 211:'pi+', 111:'pi0', 213:'rho+',
    22:'gamma'
}

def decipher(signature):
    result = ""
    reminder = signature
    if signature > 0:
        for id, name in pdgIds.items():
            while reminder % id == 0:
                if result != "":
                    result += " "
                result += name
                reminder = reminder / id
    return (result, reminder)

def print_canvas(output_name_without_extention, path, canvas=ROOT.gPad):
    if not os.path.exists(path):
        os.makedirs(path)
    canvas.Print(f"{path}/{output_name_without_extention}.png")
    canvas.Print(f"{path}/{output_name_without_extention}.pdf")
    canvas.Print(f"{path}/{output_name_without_extention}.root")

def name(pdgId):
    id = abs(pdgId)
    if id not in pdgIds:
        if id not in unknown_pdgIds:
            unknown_pdgIds.append(id)
        return "other"
    return pdgIds[id]


chain = ROOT.TChain("inclusiveMc")
for f in glob.glob("%s/*.root" % input_path):
    chain.Add(f)

all_processes = defaultdict(int)
all_matched_processes = defaultdict(int)
all_processes_bmass = defaultdict(int)
all_matched_processes_bmass = defaultdict(int)

nevents = chain.GetEntries()
nevents_bmass = 0
nevents_bmass_one_tree = 0

print(f"N selected candidates: {nevents}")

h_kpi_mass_1z  = ROOT.TH1D("h_kpi_mass_1z","Kpi mass;mass, GeV", 25, 0.816, 0.976)
h_kstar_mass_1z  = ROOT.TH1D("h_kstar_mass_1z","Kstar mass;mass, GeV", 25, 0.816, 0.976)
for event in chain:
    signature = abs(event.decay_signature)
    all_processes[signature] += 1

    if event.mc_match !=0 and event.decay_parent == event.mc_match:
        all_matched_processes[signature] += 1

    if abs(event.m-5.28) < 0.1:
        nevents_bmass += 1
        all_processes_bmass[signature] += 1

        if event.mc_match !=0 and event.decay_parent == event.mc_match:
            if signature == 443*313:
                h_kstar_mass_1z.Fill(event.kstar_mass)
            elif signature == 443*321*211:
                h_kpi_mass_1z.Fill(event.kstar_mass)
            nevents_bmass_one_tree += 1
            all_matched_processes_bmass[signature] += 1

print("All processes")
for signature, count in sorted(all_processes.items(), key=lambda x: x[1], reverse=True):
    if count > 10:
        decay, reminder = decipher(signature)
        if reminder != 1:
            decay += f" {reminder}"
        print(f"{decay:30s} {count} ({100.0*count/nevents:0.1f}%)")

print("\nProcesses where all particles originate from one decay tree")
for signature, count in sorted(all_matched_processes.items(), key=lambda x: x[1], reverse=True):
    if count > 10:
        decay, reminder = decipher(signature)
        if reminder != 1:
            decay += f" {reminder}"
        print(f"{decay:30s} {count} ({100.0*count/nevents:0.1f}%)")

print(f"\nAll processes with abs(m-5.28)<0.1: {nevents_bmass}")
shown_fraction = 0
for signature, count in sorted(all_processes_bmass.items(), key=lambda x: x[1], reverse=True):
    fraction = 100.0*count/nevents_bmass
    if fraction > 1.0:
        shown_fraction += fraction
        decay, reminder = decipher(signature)
        if reminder != 1:
            decay += f" {reminder}"
        print(f"{decay:30s} {count} ({fraction:0.1f}%)")
print(f"Other {100-shown_fraction:0.1f}%")

print(f"\nProcesses where all particles originate from one decay tree and abs(m-5.28)<0.1: {nevents_bmass_one_tree}")
shown_fraction = 0
for signature, count in sorted(all_matched_processes_bmass.items(), key=lambda x: x[1], reverse=True):
    fraction = 100.0*count/nevents_bmass_one_tree
    if fraction > 1.0:
        shown_fraction += fraction
        decay, reminder = decipher(signature)
        if reminder != 1:
            decay += f" {reminder}"
        print(f"{decay:30s} {count} ({fraction:0.1f}%)")
print(f"Other {100-shown_fraction:0.1f}%")

# Make plots
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

h_kstar_mass_1z.Draw()
h_kstar_mass_1z.SetMinimum(0)
print(f"N events (h_kstar_mass_1z): {h_kstar_mass_1z.Integral():0.1f}")
print_canvas("h_minbias_kstar_mass_1z", output_path)

h_kpi_mass_1z.Draw()
h_kpi_mass_1z.SetMinimum(0)
print(f"N events (h_kpi_mass_1z): {h_kpi_mass_1z.Integral():0.1f}")
print_canvas("h_minbias_kpi_mass_1z", output_path)


