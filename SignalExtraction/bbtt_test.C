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
//#include "../Utils/CMS_lumi_v2.h"

#endif

TTree * load(std::string iName) { 
  TFile *lFile = new TFile(iName.c_str());
  TTree *lTree  = (TTree*) lFile->FindObjectAny("Events");
  return lTree;
}

float maximum(TH1F* h, bool LOG=false){
  if(LOG){
    if(h->GetMaximum()>1000){ return 1000.*TMath::Nint(30*h->GetMaximum()/1000.); }
    if(h->GetMaximum()>  10){ return   10.*TMath::Nint(30*h->GetMaximum()/  10.); }
    return 50*h->GetMaximum(); 
  }
  else{
    if(h->GetMaximum()>  12){ return 10.*TMath::Nint((1.20*h->GetMaximum()/10.)); }
    if(h->GetMaximum()> 1.2){ return TMath::Nint((1.50*h->GetMaximum())); }
    //return 1.6*h->GetMaximum(); 
    return 1.7*h->GetMaximum(); 
  }
}

void blind(TH1F* iH,double low,double high) {
  for(int i = 0; i < iH->GetNbinsX(); i++) {
    double mass = iH->GetBinCenter(i);
    if(mass <= high && mass >= low)
      {
	iH->SetBinContent(i,0);
	iH->SetBinError(i,0);
      }
  }
}

void bbtt_test(std::string var="mt2pileup",int nbins=8, double xmin=0, double xmax=200,std::string xtitle="mt2pileup", std::string ytitle="Events", double sigscale=1,int hist=1)
{
  double massLEdges[14]    = {-0.5,-0.45,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.0,0.05,0.1,0.2};

  //SetStyle(); gStyle->SetLineStyleString(11,"20 10");
  //setTDRStyle();
  TH1::SetDefaultSumw2(1);

  std::string dir = "/afs/cern.ch/work/j/jlawhorn/public/holding/";
  //std::string dir = "/afs/cern.ch/user/j/jlawhorn/delphes-dihiggs/CombSelection/";
  
  std::stringstream scale; scale << sigscale;
  
  //Cut definitions
  double luminosity = 3000;
  std::stringstream lumi; lumi << luminosity;
  std::string objN = "(isBBTT==1 && tauCat1==1 && tauCat2==1 && ptTau1>45 && ptTau2>45 && (bTag1==2||bTag1==3||bTag1==6||bTag1==7) && (bTag2==2||bTag2==3||bTag2==6||bTag2==7) && ptB1>30 && ptB2>30)*(abs(etaB1)<2.5 && abs(etaB2)<2.5 && abs(etaTau1)<2.1 && abs(etaTau2)<2.1 && ptTrk1>0 && ptTrk2>0)*eventWeight*"+lumi.str();
  std::string objT = "(isBBTT==1 && tauCat1==1 && tauCat2==1 && ptTau1>45 && ptTau2>45 && ptB1>30 && ptB2>30)*(abs(etaB1)<2.5 && abs(etaB2)<2.5 && abs(etaTau1)<2.1 && abs(etaTau1)<2.5 && ptTrk1>0 && ptTrk2>0)*eventWeight*"+lumi.str();
  std::string jetN = objN+"*(m_svpileup>90 && m_svpileup<120 && mBB1>90 && mBB1<140)";
  std::string jetT = objT+"*(m_svpileup>90 && m_svpileup<120 && mBB1>90 && mBB1<140)";
  //signal region
  
  //--------------------------------------------------------------------------
  
  //Get the trees
  //TTree *tt_nom = load(dir+"tt-4p-0-600-v1510_14TEV_nov.root");
  TTree *tt_nom = load(dir+"tt-BLUE.root");
  TTree *tt_test = load(dir+"tt_test.root");
  TTree *tt_wpuid = load(dir+"tt_wpuid.root");
  //TTree *tt_nom = load(dir+"tt-4p-2500-100000-v1510_14TEV.root");
  //TTree *tt_test = load(dir+"tt-4p-2500-100000-v1510_14TEV_test.root");
  //TTree *tt_wpuid = load(dir+"tt-4p-2500-100000-v1510_14TEV_puid.root");
  
  //-------------------------------------------------------------------------
  
  //Get histograms
  TCanvas *canv0 = MakeCanvas("canv", "histograms", 600, 600);
  canv0->cd();
  std::string vardraw;
  TH1F *h_nom;
  if(hist)
    h_nom = new TH1F("TTnom","",nbins,xmin,xmax);
  else 
    h_nom = new TH1F("TTnom","",13, massLEdges);
  vardraw = var+">>"+"TTnom";
  tt_nom->Draw(vardraw.c_str(),jetN.c_str());
  InitSignal(h_nom);
  TH1F *h_test;
  if(hist)
    h_test = new TH1F("TTtest","",nbins,xmin,xmax);
  else 
    h_test = new TH1F("TTtest","",13, massLEdges);
  vardraw = var+">>"+"TTtest";
  tt_test->Draw(vardraw.c_str(),objT.c_str());
  InitSignal(h_test);
  TH1F *h_puid;
  if(hist)
    h_puid = new TH1F("TTpuid","",nbins,xmin,xmax);
  else 
    h_puid = new TH1F("TTpuid","",13, massLEdges);
  vardraw = var+">>"+"TTpuid";
  tt_wpuid->Draw(vardraw.c_str(),objT.c_str());
  InitSignal(h_puid);

  delete canv0;

  Double_t error=999;
  //cout << objN << endl;
  //cout << objT << endl;
  cout << jetN << endl;
  cout << jetT << endl;
  cout << "Nominal   "  << h_nom->IntegralAndError(0,h_nom->GetNbinsX(),error) << "+/-";
  cout << error << endl; error=999; cout << h_nom->GetEntries() << endl;
  cout << "Testing   "  << h_test->IntegralAndError(0,h_test->GetNbinsX(),error) << "+/-";
  cout << error << endl; error=999; cout << h_test->GetEntries() << endl;
  cout << "Test+PU   "  << h_puid->IntegralAndError(0,h_puid->GetNbinsX(),error) << "+/-";
  cout << error << endl; error=999; cout << h_puid->GetEntries() << endl;

  Float_t sc_test=h_test->Integral();
  Float_t sc_nom=h_nom->Integral();
  Float_t sc_puid=h_puid->Integral();
  h_test->Scale(sc_nom/sc_test);
  h_puid->Scale(sc_nom/sc_puid);

  cout << "----" << endl;
  cout << "Nominal   "  << h_nom->IntegralAndError(0,h_nom->GetNbinsX(),error) << "+/-";
  cout << error << endl; error=999; cout << h_nom->GetEntries() << endl;
  cout << "Testing   "  << h_test->IntegralAndError(0,h_test->GetNbinsX(),error) << "+/-";
  cout << error << endl; error=999; cout << h_test->GetEntries() << endl;
  cout << "Test+PU   "  << h_puid->IntegralAndError(0,h_puid->GetNbinsX(),error) << "+/-";
  cout << error << endl; error=999; cout << h_puid->GetEntries() << endl;

  TCanvas *canv = MakeCanvas("canv", "histograms", 800, 600);
  canv->cd();
  
  h_test->SetLineColor(kBlue);
  h_test->SetLineWidth(2.0);
  h_puid->SetLineColor(kRed);
  h_puid->SetLineWidth(2.0);
  h_nom->SetLineColor(kBlack);
  h_nom->SetLineWidth(2.0);

  h_nom->SetMaximum(std::max(maximum(h_puid,0), maximum(h_nom,0)));
  h_nom->GetXaxis()->SetTitle(xtitle.c_str());
  h_nom->GetYaxis()->SetTitle(ytitle.c_str());
  h_nom->Draw("hist");
  h_test->Draw("histsame");
  h_puid->Draw("histsame");

  TLegend* leg = new TLegend(0.65, 0.65, 0.95, 0.90);
  SetLegendStyle(leg);
  leg->AddEntry(h_nom, "Nominal Selection", "L");
  leg->AddEntry(h_test, "Generator Matching", "L");
  leg->AddEntry(h_puid, "Generator w/PU-ID", "L");
  leg->Draw();

  canv->Print((var+"_test.png").c_str());

}
