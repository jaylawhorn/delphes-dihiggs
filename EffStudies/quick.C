#if !defined(__CINT__) || defined(__MAKECINT__)
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <sstream>
#include <TH1F.h>
#include <TFile.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <Rtypes.h>
#include "TGraph.h"
#include <TMath.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TAttLine.h>
#include <TPaveText.h>
#include <TColor.h>

#include "TTree.h"
#include "TCanvas.h"

#include "../Utils/HttStyles.h"
#include "../Utils/CMS_lumi_v2.h"

#endif

void quick() {

  //setTDRStyle();

  TCanvas *canv = MakeCanvas("canv", "histograms", 800, 600);
  canv->cd();

  double eff_points[10] = { 0.55, 0.57, 0.59, 0.61, 0.63, 0.65, 0.67, 0.69, 0.71, 0.73 };
  double eff_points_2[10] = { 0.55, 0.57, 0.59, 0.61, 0.63, 0.65, 0.67, 0.69, 0.71, 0.73 };

  double sig_up[10] = { 1.22484, 1.20737, 1.19085, 1.17521, 1.16035, 1.14622, 1.13276, 1.11993, 1.10767, 1.09593 };


  TF1 *fxn = new TF1("fxn", "[0]*TMath::Sqrt([1]*x+[2])/([3]*x+[4])", 0.4, 0.9);
  fxn->SetParameter(0,114.622/10.7759);
  fxn->SetParameter(1,11548);
  fxn->SetParameter(2,658);
  fxn->SetParameter(3,12.9);
  fxn->SetParameter(4,0);

  for (Int_t i=0; i<10; i++) {
    sig_up[i]*=100;
    eff_points_2[i]/=0.65;
  }

  TGraph *gr = new TGraph(10, eff_points, sig_up);
  //TGraphAsymmErrors gr2(10, eff_points, sig_cent, 0, 0, sig_down, sig_up);

  gr->SetLineColor(kRed);
  //gr->SetLineWidth(2);
  fxn->SetLineColor(kBlue);
  //fxn->SetLineWidth(2);
  //gr2->SetFillColor(kGreen);

  gr->GetXaxis()->SetTitle("Hadronic #tau Efficiency");
  gr->GetYaxis()->SetTitle("Uncertainty on #sigma_{HH}");
  gr->Draw("alp");


  cout << fxn->Eval(0.65) << endl;

  fxn->Draw("same R");


}
