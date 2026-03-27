#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TAxis.h"

#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooJohnson.h"
#include "RooPlot.h"
#include "RooFit.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooFitResult.h"
#include "RooHist.h"

using namespace RooFit;

void RooJohnson_toyStudy() {
  gStyle->SetOptStat(0);

  // Observable and true parameters
  RooRealVar mass("mass","mass", 2.8, 3.4);

  const double trueMu    = 3.1;
  const double trueSigma = 0.05;

  RooRealVar mu_true   ("mu_true",   "mu_true",   trueMu);
  RooRealVar sigma_true("sigma_true","sigma_true",trueSigma);
  mu_true.setConstant(true);
  sigma_true.setConstant(true);

  // True Gaussian resolution model
  RooGaussian gaus_true("gaus_true","true Gaussian", mass, mu_true, sigma_true);

  // Generate a large toy sample from a pure Gaussian
  const int nEvents = 1000000;
  RooDataSet* data = gaus_true.generate(mass, nEvents);

  // --------------------------------------------------
  // Fit with a Gaussian
  // --------------------------------------------------
  RooRealVar muG   ("muG",   "muG",   3.1, 3.0, 3.2);
  RooRealVar sigmaG("sigmaG","sigmaG",0.05,0.01,0.08);

  RooGaussian gaus_fit("gaus_fit","Gaussian fit", mass, muG, sigmaG);

  std::unique_ptr<RooFitResult> rG(
    gaus_fit.fitTo(*data, Save(true), PrintLevel(-1))
  );

  // --------------------------------------------------
  // Fit with a Johnson S_U
  // --------------------------------------------------
  RooRealVar muJ    ("muJ",    "muJ",    3.1, 3.0, 3.2);
  RooRealVar lambdaJ("lambdaJ","lambdaJ",0.05, 0.001, 10.0); // >0
  RooRealVar gammaJ ("gammaJ", "gammaJ", 0.0, -1.0, 1.0);
  RooRealVar deltaJ ("deltaJ", "deltaJ", 1.0, 0.1, 100.0);   // >0

  RooJohnson johnson_fit("johnson_fit","Johnson S_U fit",
                         mass, muJ, lambdaJ, gammaJ, deltaJ);

  std::unique_ptr<RooFitResult> rJ(
    johnson_fit.fitTo(*data, Save(true), PrintLevel(-1))
  );

  // --------------------------------------------------
  // Compute chi2 for each model
  // --------------------------------------------------
  RooPlot* frameG = mass.frame(Title("Gaussian toys, Gaussian fit"));
  data->plotOn(frameG, Name("dataG"));
  gaus_fit.plotOn(frameG, LineColor(kBlue), LineWidth(2), Name("gaus_curve"));
  int nparG = rG->floatParsFinal().getSize();
  double chi2G = frameG->chiSquare("gaus_curve","dataG", nparG);

  RooPlot* frameJ = mass.frame(Title("Gaussian toys, Johnson S_{U} fit"));
  data->plotOn(frameJ, Name("dataJ"));
  johnson_fit.plotOn(frameJ, LineColor(kRed), LineWidth(2), Name("johnson_curve"));
  int nparJ = rJ->floatParsFinal().getSize();
  double chi2J = frameJ->chiSquare("johnson_curve","dataJ", nparJ);

  // Print fit summaries
  std::cout << "\n=== Gaussian fit to Gaussian toys ===\n";
  rG->Print("v");
  std::cout << "chi2/ndf (Gaussian) = " << chi2G << "\n";

  std::cout << "\n=== Johnson S_U fit to Gaussian toys ===\n";
  rJ->Print("v");
  std::cout << "chi2/ndf (Johnson) = " << chi2J << "\n\n";

  // --------------------------------------------------
  // Visual comparison on one canvas
  // --------------------------------------------------
  TCanvas* c = new TCanvas("c","Gaussian toys: Gaussian vs Johnson fit", 900, 700);
  c->Divide(1,2);

  // Top: data with both fits overlaid
  c->cd(1);
  RooPlot* frame = mass.frame(Title("Gaussian toys: Gaussian and Johnson fits"));
  data->plotOn(frame, Name("data"));
  gaus_fit.plotOn(frame,
                  LineColor(kBlue),
                  LineWidth(2),
                  Name("gaus_curve_all"));
  johnson_fit.plotOn(frame,
                     LineColor(kRed),
                     LineStyle(kDashed),
                     LineWidth(2),
                     Name("johnson_curve_all"));

  frame->GetYaxis()->SetTitle("Events / bin");
  frame->Draw();

  {
    TLegend* leg = new TLegend(0.55,0.65,0.88,0.88);
    leg->SetBorderSize(0);
    leg->AddEntry(frame->findObject("data"),               "Gaussian toys", "p");
    leg->AddEntry(frame->findObject("gaus_curve_all"),
                  Form("Gaussian fit (chi2/ndf = %.2f)", chi2G),
                  "l");
    leg->AddEntry(frame->findObject("johnson_curve_all"),
                  Form("Johnson fit (chi2/ndf = %.2f)", chi2J),
                  "l");
    leg->Draw();
  }

  // Bottom: pulls for Johnson fit
  c->cd(2);
  RooPlot* framePull = mass.frame(Title("Pulls for Johnson S_{U} fit"));
  data->plotOn(framePull, Name("data_for_pull"));
  johnson_fit.plotOn(framePull, Name("johnson_for_pull"));
  RooHist* hpull = framePull->pullHist("data_for_pull","johnson_for_pull");
  RooPlot* pullFrame = mass.frame();
  pullFrame->addPlotable(hpull,"P");
  pullFrame->GetYaxis()->SetTitle("Pull");
  pullFrame->Draw();

  c->Update();
}
