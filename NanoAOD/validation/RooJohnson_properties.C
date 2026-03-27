#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TAxis.h"

#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooJohnson.h"
#include "RooPlot.h"
#include "RooFit.h"

using namespace RooFit;

void RooJohnson_properties() {
  gStyle->SetOptStat(0);

  // Observable
  RooRealVar mass("mass","mass", 2.8, 3.4);

  // Parameters (Johnson)
  RooRealVar mu    ("mu",    "mu",    3.1,  3.0, 3.2);
  RooRealVar lambda("lambda","lambda",0.05, 1e-4, 0.5); // > 0
  RooRealVar gamma ("gamma", "gamma", 0.0, -6.0, 6.0);
  RooRealVar delta ("delta", "delta", 1.0, 1e-3, 10.0); // > 0

  lambda.setConstant(true);

  // Johnson S_U pdf
  RooJohnson johnson("johnson","Johnson S_U", mass, mu, lambda, gamma, delta);

  // Gaussian with same mu and similar scale
  RooRealVar sigma("sigma","sigma",0.05, 1e-4, 0.5);
  sigma.setConstant(true);
  RooGaussian gauss("gauss","Gaussian", mass, mu, sigma);

  // ------------------------------------------------------------
  // 1) RooGaussian vs RooJohnson with gamma = 0 delta = 1
  // ------------------------------------------------------------
  TCanvas* c1 = new TCanvas("c1","1) Gaussian vs Johnson (gamma=0 delta=1)", 800, 600);

  gamma.setVal(0.0);
  delta.setVal(1.0);

  RooPlot* frame1 = mass.frame(Title("Gaussian vs Johnson S_{U} (#gamma=0 #delta=1)"));

  gauss.plotOn(frame1,
               LineColor(kBlack),
               LineWidth(2),
               Name("gauss"));
  johnson.plotOn(frame1,
                 LineColor(kRed),
                 LineWidth(2),
                 LineStyle(kDashed),
                 Name("johnson_g0_d1"));

  frame1->GetYaxis()->SetTitle("Probability density");
  frame1->Draw();

  {
    TLegend* leg = new TLegend(0.65,0.70,0.88,0.88);
    leg->SetBorderSize(0);
    leg->AddEntry(frame1->findObject("gauss"),         "Gaussian", "l");
    leg->AddEntry(frame1->findObject("johnson_g0_d1"), "Johnson S_{U}", "l");
    leg->Draw();
  }

  c1->Update();

  // ------------------------------------------------------------
  // 2) Scan for gamma (keep delta = 1 Gaussian shown)
  // ------------------------------------------------------------
  TCanvas* c2 = new TCanvas("c2","2) Gamma scan", 800, 600);

  delta.setVal(1.0);

  RooPlot* frame2 = mass.frame(Title("Scan in #gamma (#mu=3.1, #lambda=0.05, #delta=1)"));

  // Gaussian reference
  gauss.plotOn(frame2,
               LineColor(kBlack),
               LineWidth(2),
               Name("gauss2"));

  double gammaVals[] = {-2, -1, 0.0, 1, 2};
  int gammaCols[]    = {kBlue+2, kBlue, kRed, kMagenta+2, kGreen+2};

  for (int i = 0; i < 5; ++i) {
    gamma.setVal(gammaVals[i]);
    TString name = Form("johnson_gamma_%+.1f", gammaVals[i]);
    int style = kSolid;
    if (i < 2 ) style = kDashed;
    if (i > 2 ) style = 9;
    
    johnson.plotOn(frame2,
                   LineColor(gammaCols[i]),
                   LineWidth(2),
                   LineStyle(style),
                   Name(name));
  }

  frame2->GetYaxis()->SetTitle("Probability density");
  frame2->Draw();

  {
    TLegend* leg = new TLegend(0.65,0.55,0.88,0.88);
    leg->SetBorderSize(0);
    leg->AddEntry(frame2->findObject("gauss2"), "Gaussian", "l");
    for (int i = 0; i < 5; ++i) {
      TString name = Form("johnson_gamma_%+.1f", gammaVals[i]);
      leg->AddEntry(frame2->findObject(name),
                    Form("Johnson, #gamma = %.1f", gammaVals[i]),
                    "l");
    }
    leg->Draw();
  }

  c2->Update();

  // ------------------------------------------------------------
  // 3) Scan for delta (keep gamma = 0 Gaussian shown)
  // ------------------------------------------------------------
  TCanvas* c3 = new TCanvas("c3","3) Delta scan", 800, 600);

  gamma.setVal(0.0);

  RooPlot* frame3 = mass.frame(Title("Scan in #delta (#mu=3.1, #lambda=0.05, #gamma=0)"));

  // Gaussian reference
  gauss.plotOn(frame3,
               LineColor(kBlack),
               LineWidth(2),
               Name("gauss3"));

  double deltaVals[] = {0.5, 1.0, 1.5, 2.0};
  int deltaCols[]    = {kBlue, kRed, kMagenta+2, kGreen+2};

  for (int i = 0; i < 4; ++i) {
    delta.setVal(deltaVals[i]);
    TString name = Form("johnson_delta_%.1f", deltaVals[i]);
    johnson.plotOn(frame3,
                   LineColor(deltaCols[i]),
                   LineWidth(2),
                   LineStyle(i == 1 ? kSolid : kDashed),
                   Name(name));
  }

  frame3->GetYaxis()->SetTitle("Probability density");
  frame3->Draw();

  {
    TLegend* leg = new TLegend(0.65,0.55,0.88,0.88);
    leg->SetBorderSize(0);
    leg->AddEntry(frame3->findObject("gauss3"), "Gaussian", "l");
    for (int i = 0; i < 4; ++i) {
      TString name = Form("johnson_delta_%.1f", deltaVals[i]);
      leg->AddEntry(frame3->findObject(name),
                    Form("Johnson, #delta = %.1f", deltaVals[i]),
                    "l");
    }
    leg->Draw();
  }

  c3->Update();
}
