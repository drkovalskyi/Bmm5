#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TF1.h"

void fit_efficiency()
{
  TFile *_file0 = TFile::Open("/afs/cern.ch/user/d/dmytro/www/public_html/plots/Bmm/AN/trigger_info_tmp//h_eff_eta1.0_Bsmm-2018_bsmm_HLT_DoubleMu4_3_Bs.root");
  TCanvas* c1 = (TCanvas*)_file0->Get("c1");
  c1->Draw();
  TH1* h = (TH1*)c1->GetListOfPrimitives()->FindObject("h_eff_eta1.0_Bsmm-2018_bsmm_HLT_DoubleMu4_3_Bs");
  TF1* merf = new TF1("merf", "[2]*(TMath::Erf((x - [0])*[1])+1)", 4 , 20);
  merf->SetParLimits(0,4,10);
  merf->SetParLimits(1,0,1);

  TF1* merf2 = new TF1("merf2", "[2]*([3]*(TMath::Erf((x - [0])*[1])+1) + (1-[3])*(TMath::Erf((x - [4])*[5])+1))", 4 , 20);
  merf2->SetParLimits(0,4,6);
  merf2->SetParLimits(1,0,2);
  merf2->SetParLimits(3,0,1);
  merf2->SetParLimits(4,4,10);
  merf2->SetParLimits(5,0,1);

  
  h->Fit("merf2");
  h->Draw();
}
