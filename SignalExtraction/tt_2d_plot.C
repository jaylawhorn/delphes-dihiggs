#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <sstream>
#include <TH1F.h>
#include <TH2F.h>
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

void tt_2d_plot()
{
  //variables to draw
  std::string var1= "mt2pileup:ptBB1";

  //plotting style
  SetStyle(); gStyle->SetLineStyleString(11,"20 10");
  gStyle->SetPalette(1,0);
  //histogram info
  int nbins = 24;
  double xmin = 0;
  double xmax = 600;
  int ynbins = 24;
  double ymin = 0;
  double ymax = 600;
  std::string xtitle = "m_{T2}";
  std::string ytitle = "p_{T}(bb)";
  //cuts
  double luminosity = 3000;
  std::stringstream lumi; lumi << luminosity;
  std::string objcut = "(tauCat1==1 && tauCat2==1 && ptTau1>45 && ptTau2>45 && ptB1>30 && ptB2>30 && (bTag1==2||bTag1==3||bTag1==6||bTag1==7) && (bTag2==2||bTag2==3||bTag2==6||bTag2==7))*(abs(etaTau1)<2.1 && abs(etaTau2)<2.1 && abs(etaB1)<2.5 && abs(etaB2)<2.5)*(ptTrk1>0 && ptTrk2>0)*(m_svpileup>70 && m_svpileup<120)*(mBB1>80&&mBB1<140)*(mHH>300)";
  std::string jetcut = objcut+"*eventWeight*"+lumi.str();

  //Get the tree
  std::string dir = "/afs/cern.ch/work/j/jlawhorn/public/ntuples/";
  TTree *vbfhiggs = load(dir+"HHToTTBB_14TeV.root");
  //Get histograms

  TCanvas *canv0 = MakeCanvas("canv", "histograms", 600, 600);
  canv0->cd();
  std::string vardraw;
  TH2F *jet1 = new TH2F("jet1","",nbins,xmin,xmax,ynbins,ymin,ymax);
  vardraw = var1+">>"+"jet1";
  vbfhiggs->Draw(vardraw.c_str(),jetcut.c_str());
  jet1->GetXaxis()->SetTitle(xtitle.c_str());
  jet1->GetYaxis()->SetTitle(ytitle.c_str());
  //InitStandardHist(jet1,xtitle.c_str(),ytitle.c_str(),TColor::GetColor(222,90,106));
  //TH1F *jet2 = new TH1F("jet2","",nbins,xmin,xmax);
  //vardraw = var2+">>"+"jet2";
  //vbfhiggs2->Draw(vardraw.c_str(),vbfcut2.c_str());
  //InitStandardHist(jet2,xtitle.c_str(),ytitle.c_str(),kBlue);
  //jet1->Scale(1.0/jet1->Integral());
  //jet2->Scale(1.0/jet2->Integral());
  
  TCanvas *canv = MakeCanvas("canv", "histograms", 600, 600);
  canv->cd();
  //jet2->SetMaximum(1.5*maximum(jet2,0));
  jet1->Draw("COLZ");
  //jet1->Draw("histsame");
  //canv->RedrawAxis();
  //Adding a legend
  //TLegend* leg = new TLegend(0.53, 0.65, 0.95, 0.90);
  //SetLegendStyle(leg);
  //leg->AddEntry(jet1  , "Leading jet #eta"           , "L" );
  //leg->AddEntry(jet2, "Second leading jet #eta"      , "L" );
  //leg->AddEntry(jet1  , "Phase I, PU 140, Aged"           , "L" );
  //leg->AddEntry(jet2, "Phase I, PU 50"      , "L" );
  //leg->Draw();

  //CMS preliminary 
   // const char* dataset = "CMS Simulation,  H#rightarrow#tau#tau, 3000 fb^{-1} at 14 TeV";
   // const char* category = "";
   // CMSPrelim(dataset, "", 0.16, 0.835);
   // TPaveText* chan     = new TPaveText(0.52, 0.35, 0.91, 0.55, "tlbrNDC");
   // chan->SetBorderSize(   0 );
   // chan->SetFillStyle(    0 );
   // chan->SetTextAlign(   12 );
   // chan->SetTextSize ( 0.05 );
   // chan->SetTextColor(    1 );
   // chan->SetTextFont (   62 );
   // chan->AddText(category);
   // chan->Draw();
   //-------------------------------------------------------------------------
   //Save histograms
   canv->Print((var1+".png").c_str());
}
