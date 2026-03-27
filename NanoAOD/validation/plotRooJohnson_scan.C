#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "RooRealVar.h"
#include "RooJohnson.h"
#include "RooPlot.h"
#include "RooFit.h"
#include "TAxis.h"

using namespace RooFit;

void plotRooJohnson_scan() {
  // Observable
  RooRealVar mass("mass","mass", 2.8, 3.4);

  // Parameters
  RooRealVar mu("mu","mu", 3.1);
  RooRealVar lambda("lambda","lambda", 0.05, 1e-3, 10.0);
  RooRealVar gamma("gamma","gamma", 0.0);
  RooRealVar delta ("delta", "delta",  1.0, 1e-3, 10.0);
  
  // Johnson SU pdf
  RooJohnson johnson("johnson","Johnson S_U", mass, mu, lambda, gamma, delta);

  TCanvas* c1 = new TCanvas("c1","RooJohnson parameter scan",800,600);
  gStyle->SetOptStat(0);

  // Example 1: scan gamma
  RooPlot* frame1 = mass.frame(Title("Johnson S_{U}: scan in #gamma (#mu=3.1, #sigma=0.05)"));

  double gammaVals[] = {-2.0, -1.0, 0.0, 1.0, 2.0};
  int    gammaCols[] = {kBlue+2, kBlue, kBlack, kRed, kRed+2};

  for (int i = 0; i < 5; ++i) {
    gamma.setVal(gammaVals[i]);
    TString name = Form("g_%+.0f", gammaVals[i]);
    johnson.plotOn(frame1,
                   LineColor(gammaCols[i]),
                   LineStyle(i == 2 ? kSolid : kDashed),
                   Name(name));
  }

  frame1->GetYaxis()->SetTitle("Probability density");
  frame1->Draw();

  auto leg1 = new TLegend(0.65,0.65,0.88,0.88);
  for (int i = 0; i < 5; ++i) {
    TString name = Form("g_%+.0f", gammaVals[i]);
    leg1->AddEntry(frame1->findObject(name), Form("#gamma = %.1f", gammaVals[i]), "l");
  }
  leg1->SetBorderSize(0);
  leg1->Draw();

  c1->SetGrid();

  // Example 2: new canvas to compare different deltas at fixed gamma
  TCanvas* c2 = new TCanvas("c2","RooJohnson delta scan",800,600);

  gamma.setVal(1.0);  // fix gamma for this scan

  RooPlot* frame2 = mass.frame(Title("Johnson S_{U}: scan in #delta (#mu=3.1, #sigma=0.05)"));

  double deltaVals[] = {0.5, 1.0, 2.0, 3.0};
  int    deltaCols[] = {kBlue, kBlack, kRed, kGreen+2};

  for (int i = 0; i < 4; ++i) {
    delta.setVal(deltaVals[i]);
    TString name = Form("d_%.1f", deltaVals[i]);
    johnson.plotOn(frame2,
                   LineColor(deltaCols[i]),
                   LineStyle(i == 1 ? kSolid : kDashed),
                   Name(name));
  }

  frame2->GetYaxis()->SetTitle("Probability density");
  frame2->Draw();

  auto leg2 = new TLegend(0.65,0.65,0.88,0.88);
  for (int i = 0; i < 4; ++i) {
    TString name = Form("d_%.1f", deltaVals[i]);
    leg2->AddEntry(frame2->findObject(name), Form("#delta = %.1f", deltaVals[i]), "l");
  }
  leg2->SetBorderSize(0);
  leg2->Draw();

  c2->SetGrid();
}
