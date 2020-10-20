#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"
void plot_sensitivity()
{
  TCanvas *c1 = new TCanvas("c1", "", 500, 500);
  const Int_t n = 14;
  Double_t mva[n]               = {0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 0.992, 0.994, 0.996, 0.998};
  Double_t bsmm_significance[n] = {12.0, 12.6, 12.7, 13.3, 13.9, 14.5, 15.3, 16.2, 17.2, 18.5,  18.6,  18.4,  17.2,  14.4};
  Double_t bmm_significance[n]  = { 1.4,  1.3,  1.4,  1.4,  1.4,  1.5,  1.7,  1.8,  1.9,  2.0,   2.0,   2.0,   2.1,   1.7};
  TGraph* gr = new TGraph(n, mva, bsmm_significance);
  gr->SetMarkerStyle(20);
  gr->SetMarkerColor(kBlue);
  gr->GetXaxis()->SetTitle("MVA");
  gr->GetYaxis()->SetTitle("#sigma(Bs)");
  gr->Draw("ACP");
  
}
