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


//#include "HttStyles.h"
#include "HttStyles.cc"

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


void test_msv(std::string var="m_sv",int nbins=25, double xmin=0, double xmax=250,std::string xtitle="m_sv", std::string ytitle="Events")
{
  SetStyle(); gStyle->SetLineStyleString(11,"20 10");
  TH1::SetDefaultSumw2(1);

  std::string dir = "/afs/cern.ch/work/a/arapyan/public/svfitsamples/";
  double sigscale = 10;
  double sigscale1 = 10; 
  std::stringstream scale; scale << sigscale;
  std::stringstream scale1; scale1 << sigscale1;

  //Cut definitions
  double luminosity = 3000;
  std::stringstream lumi; lumi << luminosity;
  std::string objcut = "(ptTau1>30 && ptTau2>30 && abs(etaTau1) <4.0 && abs(etaTau1)<4.0 )";
  //std::string jetcut = objcut+"*(ptJet1>30 && ptJet2>30 && abs(etaJet1) <4.7 && abs(etaJet2) <4.7 )";
  std::string mthcut = objcut+"*( ( (tauCat1==3 && tauCat2==2) || (tauCat2==3 && tauCat1==2) ) )";

  //signal region
  //std::string vbfcut = ththcut+"*eventWeight*(eventType==0)*"+lumi.str();
  //std::string zttcut = ththcut+"*eventWeight*(eventType==1)*"+lumi.str();
  //std::string ttbarcut = ththcut+"*eventWeight*(eventType==3)*"+lumi.str();
  //std::string ewkcut = ththcut+"*eventWeight*(eventType==2)*"+lumi.str();
  //std::string othercut = ththcut+"*eventWeight*(eventType==4 || eventType==2)*"+lumi.str();
  
  //--------------------------------------------------------------------------
  
  //Get the trees
  TTree *tree = load(dir+"hh.root"); 
  TTree *tree2 = load(dir+"Bjj-vbf.root"); 
  TTree *tree3 = load(dir+"vbf_bgd.root"); 

  //-------------------------------------------------------------------------
  
  //Get histograms
  TCanvas *canv0 = MakeCanvas("canv", "histograms", 600, 600);
  canv0->cd();
  std::string vardraw;
  TH1F *vbf = new TH1F("VBFH","",nbins,xmin,xmax);
  vardraw = var+">>"+"VBFH";
  tree->Draw(vardraw.c_str(),mthcut.c_str());
  InitSignal(vbf);
  vbf->SetLineColor(kBlack);
  //TH1F *ttbar = new TH1F("TTbar","",nbins,xmin,xmax);
  //vardraw = var+">>"+"TTbar";
  //tree->Draw(vardraw.c_str(),ttbarcut.c_str());
  //InitHist(ttbar, xtitle.c_str(), ytitle.c_str(), TColor::GetColor(155,152,204), 1001);
  TH1F *ztt = new TH1F("Ztt","",nbins,xmin,xmax);
  vardraw = var+">>"+"Ztt";
  tree2->Draw(vardraw.c_str(),mthcut.c_str());
  InitHist(ztt, xtitle.c_str(), ytitle.c_str(), TColor::GetColor(248,206,104), 1001);
  //TH1F *ewk = new TH1F("Ewk","",nbins,xmin,xmax);
  //vardraw = var+">>"+"Ewk";
  //tree->Draw(vardraw.c_str(),ewkcut.c_str());
  //InitHist(ewk, xtitle.c_str(), ytitle.c_str(),  TColor::GetColor(222,90,106), 1001);
  TH1F *other = new TH1F("Other","",nbins,xmin,xmax);
  vardraw = var+">>"+"Other";
  tree3->Draw(vardraw.c_str(),mthcut.c_str());
  InitHist(other, xtitle.c_str(), ytitle.c_str(),  TColor::GetColor(222,90,106), 1001);
  
  delete canv0;

  //----------------------------------------------------------------------------
  //Print out the yields
  /*  Double_t error=0.0;
  ofstream outfile;
  outfile.open("yields_test.txt");
  outfile << "Yields for the signal region." << std::endl;
  outfile << "VBF   "  << vbf->IntegralAndError(0,vbf->GetNbinsX(),error) << "+/-" << error << endl;
  outfile << "TTbar   "  << ttbar->IntegralAndError(0,ttbar->GetNbinsX(),error) << "+/-" << error << endl;
  outfile << "Ztt    "  << ztt->IntegralAndError(0,ztt->GetNbinsX(),error) << "+/-" << error << endl;
  //outfile << "ewk    "  << ewk->IntegralAndError(0,ewk->GetNbinsX(),error) << "+/-" << error << endl;
  outfile << "other   "  << other->IntegralAndError(0,other->GetNbinsX(),error) << "+/-" << error << endl;
  outfile << "S/sqrt(B)    "  << vbf->Integral()/(other->Integral()+ztt->Integral()+ttbar->Integral()) << endl;
  //--------------------------------------------------------------------------
  //continue outputing
  //outfile << "Ewk total    "  << ewk->IntegralAndError(0,ewk->GetNbinsX(),error) << "+/-" << error << endl;
  outfile << endl << endl << endl;
  outfile << "In the signal region (100,150GeV)  " <<endl;
  outfile << "VBF   "  << vbf->IntegralAndError(5,11,error) << "+/-" << error << endl;
  outfile << "TTbar    "  << ttbar->IntegralAndError(5,11,error) << "+/-" << error << endl;
  outfile << "Ztt    "  << ztt->IntegralAndError(5,11,error) << "+/-" << error << endl;
  //outfile << "ewk    "  << ewk->IntegralAndError(5,11,error) << "+/-" << error << endl;
  outfile << "other    "  << other->IntegralAndError(5,11,error) << "+/-" << error << endl;
  outfile.close();*/
  //-----------------------------------------------------------------------
  //Draw the histograms
  TCanvas *canv = MakeCanvas("canv", "histograms", 600, 600);
  canv->cd();
  //ewk->Add(other); ttbar->Add(ewk); 
  //ztt->Add(ttbar); vbf->Add(ztt);
  other->Add(ztt);
  //Error band stat
  TH1F* errorBand = (TH1F*)vbf ->Clone("errorBand");
  errorBand  ->SetMarkerSize(0);
  errorBand  ->SetFillColor(13);
  errorBand  ->SetFillStyle(3013);
  errorBand  ->SetLineWidth(1);
  //  for(int idx=0; idx<errorBand->GetNbinsX(); ++idx){
  //     if(errorBand->GetBinContent(idx)>0){
  //       std::cout << "Uncertainties on summed background samples: " << errorBand->GetBinError(idx)/errorBand->GetBinContent(idx) << std::endl;
  //       break;
  //     }
  //}
  Float_t n=other->GetEntries();
  other->Scale(1.0/n);
  n=vbf->GetEntries();
  vbf->Scale(1.0/n);
  //for (Int_t i=1; i<nbins+1; i++) {
  //cout << "i: " << i << " " <<  other->GetBinCenter(i) << " " << other->Integral(0,i) << " " << other->Integral(i,nbins+1) << endl;
  //}
  cout << " other " << other->Integral(10,15) << endl;
  cout << " htt " << vbf->Integral(10,15) << endl;
  //n=ztt->GetEntries();
  //ztt->Scale(1.0/n);
  other->SetMaximum(1.2*std::max(maximum(vbf, 0), maximum(other, 0)));
  //blind(vbf,75,150);
  //ttbar->Draw("histsame");
  //ewk->Draw("histsame");
  other->Draw("hist");
  //ztt->Draw("histsame");
  vbf->Draw("histsame");
  //errorBand->Draw("e2same");
  canv->RedrawAxis();
  //---------------------------------------------------------------------------
  //Adding a legend
  TLegend* leg = new TLegend(0.53, 0.65, 0.95, 0.90);
  SetLegendStyle(leg);
  leg->AddEntry(vbf , "HH"             , "F");
  //leg->AddEntry(ztt  , "bjj-vbf"  , "F" );
  //leg->AddEntry(ttbar, "t#bar{t}"              , "F" );
  //leg->AddEntry(ewk  , "Electroweak"           , "F" );
  leg->AddEntry(other, "Other"                 , "F" );
  //leg->AddEntry(errorBand,"bkg. uncertainty","F");
  leg->Draw();
  //---------------------------------------------------------------------------

  //CMS preliminary 
  const char* dataset = "CMS Preliminary, H#rightarrow#tau#tau, 3.0 ab^{-1} at 14 TeV";
  const char* category = "";
  CMSPrelim(dataset, "#mu#tau_{h}", 0.17, 0.835);
  //CMSPrelim(dataset, "", 0.16, 0.835);
  TPaveText* chan     = new TPaveText(0.52, 0.35, 0.91, 0.55, "tlbrNDC");
  chan->SetBorderSize(   0 );
  chan->SetFillStyle(    0 );
  chan->SetTextAlign(   12 );
  chan->SetTextSize ( 0.05 );
  chan->SetTextColor(    1 );
  chan->SetTextFont (   62 );
  chan->AddText(category);
  chan->Draw();
  //-------------------------------------------------------------------------
  //Save histograms
  canv->Print((var+"_test.png").c_str());

}
