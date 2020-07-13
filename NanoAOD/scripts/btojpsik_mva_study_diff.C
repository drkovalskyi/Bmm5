#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include <vector>
#include <string>
#include "TLegend.h"

using namespace std;

string path = "/afs/cern.ch/user/d/dmytro/www/public_html/plots/bmm5_NanoAODv6-507/bjpsik-fit/";

struct Plot{
  string f1;
  EColor color1;
  string label1;
  string f2;
  EColor color2;
  string label2;
  string name;
};

vector<Plot> plots {
  { "Run2016B-alpha_eff_reweighted.root", kBlue, "Run2016B", "Run2016H-alpha_eff_reweighted.root", kRed, "Run2016H", "Run2016BvsRun2016H-alpha_eff_reweighted"},
  { "Run2016B-bdt_eff_reweighted.root",   kBlue, "Run2016B", "Run2016H-bdt_eff_reweighted.root",   kRed, "Run2016H", "Run2016BvsRun2016H-bdt_eff_reweighted"},
  { "Run2016B-mva_eff_reweighted.root",   kBlue, "Run2016B", "Run2016H-mva_eff_reweighted.root",   kRed, "Run2016H", "Run2016BvsRun2016H-mva_eff_reweighted"},
  { "Run2016B-iso_eff_reweighted.root",   kBlue, "Run2016B", "Run2016H-iso_eff_reweighted.root",   kRed, "Run2016H", "Run2016BvsRun2016H-iso_eff_reweighted"},
  { "Run2016B-sl3d_eff_reweighted.root",  kBlue, "Run2016B", "Run2016H-sl3d_eff_reweighted.root",  kRed, "Run2016H", "Run2016BvsRun2016H-sl3d_eff_reweighted"},
  { "Run2016B-pvip_eff_reweighted.root",  kBlue, "Run2016B", "Run2016H-pvip_eff_reweighted.root",  kRed, "Run2016H", "Run2016BvsRun2016H-pvip_eff_reweighted"},
  { "Run2016B-m1iso_eff_reweighted.root", kBlue, "Run2016B", "Run2016H-m1iso_eff_reweighted.root", kRed, "Run2016H", "Run2016BvsRun2016H-m1iso_eff_reweighted"},
  { "Run2016B-m2iso_eff_reweighted.root", kBlue, "Run2016B", "Run2016H-m2iso_eff_reweighted.root", kRed, "Run2016H", "Run2016BvsRun2016H-m2iso_eff_reweighted"},
  { "Run2016B-spvip_eff_reweighted.root", kBlue, "Run2016B", "Run2016H-spvip_eff_reweighted.root", kRed, "Run2016H", "Run2016BvsRun2016H-spvip_eff_reweighted"},
    };

void make_plot(Plot plot)
{
  TFile* f1 = TFile::Open((path + plot.f1).c_str());
  TCanvas* c1 = (TCanvas*)f1->Get("c1");
  TGraph* gr1 = (TGraph*)c1->FindObject("Graph");
  gr1->SetMarkerColor(plot.color1);
  TFile* f2 = TFile::Open((path + plot.f2).c_str());
  TCanvas* c2 = (TCanvas*)f2->Get("c1");
  TGraph* gr2 = (TGraph*)c2->FindObject("Graph");
  gr2->SetMarkerColor(kRed);
  TCanvas cc("cc","",800,800);
  cc.SetGridx();
  cc.SetGridy();
  gr1->Draw("AP*");
  gr2->Draw("P");
  TLegend legend(0.10,0.80,0.30,0.90);
  legend.AddEntry(gr1, plot.label1.c_str());
  legend.AddEntry(gr2, plot.label2.c_str());
  legend.Draw();
  cc.Print((path + plot.name + ".png").c_str());
}

void btojpsik_mva_study_diff(){
  for (auto plot: plots){
    make_plot(plot);
  }
}
