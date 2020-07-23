#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include <vector>
#include <string>
#include "TLegend.h"

using namespace std;

string path = "/afs/cern.ch/user/d/dmytro/www/public_html/plots/bmm5_NanoAODv6-507/bjpsik-fit/";

struct Hist{
  string f;
  EColor color;
  string label;
};

struct Plot{
  vector<Hist> hists;
  string name;
};

vector<Plot> plots {
  { { {"Run2016B-alpha_eff_reweighted.root", kBlue, "Run2016B"}, {"Run2016H-alpha_eff_reweighted.root", kRed, "Run2016H"} }, "Run2016BvsRun2016H-alpha_eff_reweighted"},
  { { {"Run2016B-bdt_eff_reweighted.root",   kBlue, "Run2016B"}, {"Run2016H-bdt_eff_reweighted.root",   kRed, "Run2016H"} }, "Run2016BvsRun2016H-bdt_eff_reweighted"},
  { { {"Run2016B-mva_eff_reweighted.root",   kBlue, "Run2016B"}, {"Run2016H-mva_eff_reweighted.root",   kRed, "Run2016H"} }, "Run2016BvsRun2016H-mva_eff_reweighted"},
  { { {"Run2016B-iso_eff_reweighted.root",   kBlue, "Run2016B"}, {"Run2016H-iso_eff_reweighted.root",   kRed, "Run2016H"} }, "Run2016BvsRun2016H-iso_eff_reweighted"},
  { { {"Run2016B-sl3d_eff_reweighted.root",  kBlue, "Run2016B"}, {"Run2016H-sl3d_eff_reweighted.root",  kRed, "Run2016H"} }, "Run2016BvsRun2016H-sl3d_eff_reweighted"},
  { { {"Run2016B-pvip_eff_reweighted.root",  kBlue, "Run2016B"}, {"Run2016H-pvip_eff_reweighted.root",  kRed, "Run2016H"} }, "Run2016BvsRun2016H-pvip_eff_reweighted"},
  { { {"Run2016B-m1iso_eff_reweighted.root", kBlue, "Run2016B"}, {"Run2016H-m1iso_eff_reweighted.root", kRed, "Run2016H"} }, "Run2016BvsRun2016H-m1iso_eff_reweighted"},
  { { {"Run2016B-m2iso_eff_reweighted.root", kBlue, "Run2016B"}, {"Run2016H-m2iso_eff_reweighted.root", kRed, "Run2016H"} }, "Run2016BvsRun2016H-m2iso_eff_reweighted"},
  { { {"Run2016B-spvip_eff_reweighted.root", kBlue, "Run2016B"}, {"Run2016H-spvip_eff_reweighted.root", kRed, "Run2016H"} }, "Run2016BvsRun2016H-spvip_eff_reweighted"},

  { { {"Run2016B-alpha_eff_reweighted.root", kBlue, "Run2016B"}, {"Run2016H-alpha_eff_reweighted.root", kRed, "Run2016H"} }, "Run2016BvsRun2016H-alpha_eff_reweighted"},
  { { {"Run2016B-bdt_eff_reweighted.root",   kBlue, "Run2016B"}, {"Run2016H-bdt_eff_reweighted.root",   kRed, "Run2016H"} }, "Run2016BvsRun2016H-bdt_eff_reweighted"},
  { { {"Run2016B-mva_eff_reweighted.root",   kBlue, "Run2016B"}, {"Run2016H-mva_eff_reweighted.root",   kRed, "Run2016H"} }, "Run2016BvsRun2016H-mva_eff_reweighted"},
  { { {"Run2016B-iso_eff_reweighted.root",   kBlue, "Run2016B"}, {"Run2016H-iso_eff_reweighted.root",   kRed, "Run2016H"} }, "Run2016BvsRun2016H-iso_eff_reweighted"},
  { { {"Run2016B-sl3d_eff_reweighted.root",  kBlue, "Run2016B"}, {"Run2016H-sl3d_eff_reweighted.root",  kRed, "Run2016H"} }, "Run2016BvsRun2016H-sl3d_eff_reweighted"},
  { { {"Run2016B-pvip_eff_reweighted.root",  kBlue, "Run2016B"}, {"Run2016H-pvip_eff_reweighted.root",  kRed, "Run2016H"} }, "Run2016BvsRun2016H-pvip_eff_reweighted"},
  { { {"Run2016B-m1iso_eff_reweighted.root", kBlue, "Run2016B"}, {"Run2016H-m1iso_eff_reweighted.root", kRed, "Run2016H"} }, "Run2016BvsRun2016H-m1iso_eff_reweighted"},
  { { {"Run2016B-m2iso_eff_reweighted.root", kBlue, "Run2016B"}, {"Run2016H-m2iso_eff_reweighted.root", kRed, "Run2016H"} }, "Run2016BvsRun2016H-m2iso_eff_reweighted"},
  { { {"Run2016B-spvip_eff_reweighted.root", kBlue, "Run2016B"}, {"Run2016H-spvip_eff_reweighted.root", kRed, "Run2016H"} }, "Run2016BvsRun2016H-spvip_eff_reweighted"},

  { { {"Run2017B-alpha_eff_reweighted.root", kBlue, "Run2017B"}, {"Run2017C-alpha_eff_reweighted.root", kRed, "Run2017C"} }, "Run2017BvsRun2017C-alpha_eff_reweighted"},
  { { {"Run2017B-bdt_eff_reweighted.root",   kBlue, "Run2017B"}, {"Run2017C-bdt_eff_reweighted.root",   kRed, "Run2017C"} }, "Run2017BvsRun2017C-bdt_eff_reweighted"},
  { { {"Run2017B-mva_eff_reweighted.root",   kBlue, "Run2017B"}, {"Run2017C-mva_eff_reweighted.root",   kRed, "Run2017C"} }, "Run2017BvsRun2017C-mva_eff_reweighted"},
  { { {"Run2017B-iso_eff_reweighted.root",   kBlue, "Run2017B"}, {"Run2017C-iso_eff_reweighted.root",   kRed, "Run2017C"} }, "Run2017BvsRun2017C-iso_eff_reweighted"},
  { { {"Run2017B-sl3d_eff_reweighted.root",  kBlue, "Run2017B"}, {"Run2017C-sl3d_eff_reweighted.root",  kRed, "Run2017C"} }, "Run2017BvsRun2017C-sl3d_eff_reweighted"},
  { { {"Run2017B-pvip_eff_reweighted.root",  kBlue, "Run2017B"}, {"Run2017C-pvip_eff_reweighted.root",  kRed, "Run2017C"} }, "Run2017BvsRun2017C-pvip_eff_reweighted"},
  { { {"Run2017B-m1iso_eff_reweighted.root", kBlue, "Run2017B"}, {"Run2017C-m1iso_eff_reweighted.root", kRed, "Run2017C"} }, "Run2017BvsRun2017C-m1iso_eff_reweighted"},
  { { {"Run2017B-m2iso_eff_reweighted.root", kBlue, "Run2017B"}, {"Run2017C-m2iso_eff_reweighted.root", kRed, "Run2017C"} }, "Run2017BvsRun2017C-m2iso_eff_reweighted"},
  { { {"Run2017B-spvip_eff_reweighted.root", kBlue, "Run2017B"}, {"Run2017C-spvip_eff_reweighted.root", kRed, "Run2017C"} }, "Run2017BvsRun2017C-spvip_eff_reweighted"},

  { { {"Run2018A-alpha_eff_reweighted.root", kBlue, "Run2018A"}, {"Run2018D-alpha_eff_reweighted.root", kRed, "Run2018D"} }, "Run2018AvsRun2018D-alpha_eff_reweighted"},
  { { {"Run2018A-bdt_eff_reweighted.root",   kBlue, "Run2018A"}, {"Run2018D-bdt_eff_reweighted.root",   kRed, "Run2018D"} }, "Run2018AvsRun2018D-bdt_eff_reweighted"},
  { { {"Run2018A-mva_eff_reweighted.root",   kBlue, "Run2018A"}, {"Run2018D-mva_eff_reweighted.root",   kRed, "Run2018D"} }, "Run2018AvsRun2018D-mva_eff_reweighted"},
  { { {"Run2018A-iso_eff_reweighted.root",   kBlue, "Run2018A"}, {"Run2018D-iso_eff_reweighted.root",   kRed, "Run2018D"} }, "Run2018AvsRun2018D-iso_eff_reweighted"},
  { { {"Run2018A-sl3d_eff_reweighted.root",  kBlue, "Run2018A"}, {"Run2018D-sl3d_eff_reweighted.root",  kRed, "Run2018D"} }, "Run2018AvsRun2018D-sl3d_eff_reweighted"},
  { { {"Run2018A-pvip_eff_reweighted.root",  kBlue, "Run2018A"}, {"Run2018D-pvip_eff_reweighted.root",  kRed, "Run2018D"} }, "Run2018AvsRun2018D-pvip_eff_reweighted"},
  { { {"Run2018A-m1iso_eff_reweighted.root", kBlue, "Run2018A"}, {"Run2018D-m1iso_eff_reweighted.root", kRed, "Run2018D"} }, "Run2018AvsRun2018D-m1iso_eff_reweighted"},
  { { {"Run2018A-m2iso_eff_reweighted.root", kBlue, "Run2018A"}, {"Run2018D-m2iso_eff_reweighted.root", kRed, "Run2018D"} }, "Run2018AvsRun2018D-m2iso_eff_reweighted"},
  { { {"Run2018A-spvip_eff_reweighted.root", kBlue, "Run2018A"}, {"Run2018D-spvip_eff_reweighted.root", kRed, "Run2018D"} }, "Run2018AvsRun2018D-spvip_eff_reweighted"},
    };

void make_plot(Plot plot)
{
  vector<TGraph*> graphs;
  for (auto hist: plot.hists){
    TFile* f = TFile::Open((path + hist.f).c_str());
    TCanvas* c1 = (TCanvas*)f->Get("c1");
    TGraph* gr = (TGraph*)c1->FindObject("Graph");
    gr->SetMarkerColor(hist.color);
    graphs.push_back(gr);
  }
  TCanvas cc("cc","",800,800);
  cc.SetGridx();
  cc.SetGridy();
  for (unsigned int i=0; i < graphs.size(); ++i){
    if (i == 0)
      graphs.at(i)->Draw("AP*");
    else
      graphs.at(i)->Draw("P");
  }
  TLegend legend(0.10, 0.80, 0.30, 0.90);
  for (unsigned int i=0; i < graphs.size(); ++i){
    legend.AddEntry(graphs.at(i), plot.hists.at(i).label.c_str());
  }
  legend.Draw();
  cc.Print((path + plot.name + ".png").c_str());
}

void btojpsik_mva_study_diff(){
  for (auto plot: plots){
    make_plot(plot);
  }
}
