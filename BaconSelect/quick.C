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

//#include "../Utils/tdrstyle.h"
#include "../delphes-dihiggs/Utils/HttStyles.h"
#include "../delphes-dihiggs/Utils/CMS_lumi_v2.h"

#endif
void quick () {

  std::string moredelph="*(1)*(mTT>90&&mTT<120)*(mBB1<140&&mBB1>90)";//*(mTT>90 && mTT<120)";
  std::string morebacon="*(1)*(mTT>90&&mTT<120)*(mBB1>90 && mBB1<140)";//*(mTT>90 && mTT<120)";

  std::string var="mt2pileup";

  std::string xtitle="m_{T2}";
  std::string ytitle="Events";
  std::string filename="mT2_hh-12-02.png";

  TFile *delphes = new TFile("/data/blue/Bacon/029a/Upgrade/oct-5-bbtt/gFHHTobbtautau.root", "READ");
  TFile *bacon = new TFile("/data/blue/jlawhorn/hh-bacon-feb-11.root", "READ");

  //TFile *delphes = new TFile("/data/blue/Bacon/029a/Upgrade/sep-2-bbtt/tt.root", "READ");
  //TFile *bacon = new TFile("/data/blue/jlawhorn/tt-bacon-feb-11.root", "READ");

  TH1::SetDefaultSumw2(1);
  setTDRStyle();

  //TFile *delphes = new TFile("/data/blue/jlawhorn/fs_delphes_tt_just_b.root", "READ");
  //TFile *delphes = new TFile("/data/blue/Bacon/029a/Upgrade/sep-2-bbtt/tt-4p-0-600-v1510_14TEV.root", "READ");
  //TFile *d2 = new TFile("tt-4p-0-600-v1510_14TEV_test.root", "READ");
  //TFile *bacon = new TFile("hh_proc.root", "READ");
  //TFile *bacon = new TFile("/data/blue/jlawhorn/fs_bacon_tt_just_b.root", "READ");
  //TFile *bacon = new TFile("lep_cleaning.root", "READ");

  TTree *dTree = (TTree *) delphes->Get("Events");
  TTree *bTree = (TTree *) bacon->Get("Events");
  
  TH1D *dMtt = new TH1D("dMtt", "", 10, 0, 600);
  TH1D *bMtt = new TH1D("bMtt", "", 10, 0, 600);
  TH1D *bMtt2 = new TH1D("bMtt2", "", 10, 0, 600);
  
  std::string delphescut = "(tauCat1==1 && tauCat2==1 && ptTau1>45 && ptTau2>45 && ptB1>30 && ptB2>30 && (bTag1==2||bTag1==3||bTag1==6||bTag1==7) && (bTag2==2||bTag2==3||bTag2==6||bTag2==7))*(abs(etaTau1)<2.1 && abs(etaTau2)<2.1 && abs(etaB1)<2.5 && abs(etaB2)<2.5)*eventWeight*3000*(mBB1<200)";
  
  //std::string baconcut = "(tauCat1==1 && tauCat2==1 && (ptTau1>45 && ptTau2>45) && ptB1>30 && ptB2>30 && abs(etaTau1_gen)<2.1 && abs(etaTau2_gen)<2.1 && abs(etaB1_gen)<2.5 && abs(etaB2_gen)<2.5)*3000*eventWeight*0.95*0.95*0.95*0.95*(nBjetWithHadTau==0)*0.71*0.71";
  
  std::string baconcut = "(tauCat1==1 && tauCat2==1 && ptTau1>45 && ptTau2>45 && ptB1>30 && ptB2>30 && abs(etaTau1_gen)<2.1 && abs(etaTau2_gen)<2.1 && abs(etaB1_gen)<2.5 && abs(etaB2_gen)<2.5)*3000*eventWeight*0.95*0.95*0.95*0.95*(nBjetWithHadTau==0)"; //0.71

  std::string baconcut2 = "(tauCat1==1 && tauCat2==1 && ptTau1>45 && ptTau2>45 && ptB1>30 && ptB2>30 && abs(etaTau1_gen)<2.1 && abs(etaTau2_gen)<2.1 && abs(etaB1_gen)<2.5 && abs(etaB2_gen)<2.5)*3000*eventWeight*0.95*0.95*0.95*0.95*(nBjetWithHadTau>0)*0.58";

  std::string vardraw = var+">>dMtt";
  dTree->Draw(vardraw.c_str(), (delphescut+moredelph).c_str());

  //cout << (baconcut+morebacon).c_str() << endl;
  //cout << (baconcut2+morebacon).c_str() << endl;
  
  vardraw = var+">>bMtt";
  bTree->Draw(vardraw.c_str(), (baconcut+morebacon).c_str());
  
  vardraw = var+">>bMtt2";
  bTree->Draw(vardraw.c_str(), (baconcut2+morebacon).c_str());

  bMtt->Add(bMtt2);

  TCanvas *canv = MakeCanvas("canv", "histograms", 800, 600);
  canv->cd();

  dMtt->SetLineWidth(2);
  dMtt->SetLineColor(kRed);
  bMtt->SetLineWidth(2);
  bMtt->SetLineColor(kBlue);

  Double_t error=0;
  cout << "Delphes: " << endl;
  cout << dMtt->IntegralAndError(0,20,error) << " +/- ";
  cout << error << endl;
  cout << dMtt->GetEntries() << endl;
  error=0;
  cout << "Bacon " << endl;
  cout << bMtt->IntegralAndError(0,20,error) << " +/- ";
  cout << error << endl;
  cout << bMtt->GetEntries() << endl;

  /*  Float_t scale=dMtt->Integral(0,30);
  dMtt->Scale(1/scale);
  scale=bMtt->Integral(0,30);
  bMtt->Scale(1/scale);*/
  //scale=dMtt->Integral(0,20);
  //dMtt->Scale(1/scale);

  dMtt->GetXaxis()->SetTitle(xtitle.c_str());
  dMtt->GetYaxis()->SetTitle(ytitle.c_str());

  dMtt->GetYaxis()->SetRangeUser(0, 1.5*TMath::Max(dMtt->GetMaximum(), bMtt->GetMaximum()));
  dMtt->Draw("hist");
  bMtt->Draw("hist same");

  TH1F* errorBand2 = (TH1F*)dMtt->Clone("errorBand");
  errorBand2  ->SetMarkerSize(0);
  errorBand2  ->SetFillColor(kRed);
  errorBand2  ->SetFillStyle(3013);
  errorBand2  ->SetLineWidth(1);

  TH1F* errorBand = (TH1F*)bMtt->Clone("errorBand");
  errorBand  ->SetMarkerSize(0);
  errorBand  ->SetFillColor(kBlue);
  errorBand  ->SetFillStyle(3013);
  errorBand  ->SetLineWidth(1);

  errorBand->Draw("e2same");
  errorBand2->Draw("e2same");

  //Adding a legend
  TLegend* leg = new TLegend(0.65, 0.65, 0.95, 0.90);
  SetLegendStyle(leg);
  leg->AddEntry(dMtt  , "Delphes"   , "l" );
  leg->AddEntry(bMtt  , "Feb-15"     , "l" );
  leg->Draw();

  canv->SaveAs(filename.c_str());


}
