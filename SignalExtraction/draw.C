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


void draw(std::string var="mHH", std::string cut="(1)", int chan_cat=0, int nbins=25, double xmin=100, double xmax=1600,std::string xtitle="", std::string ytitle="Events")
{
  SetStyle(); gStyle->SetLineStyleString(11,"20 10");
  TH1::SetDefaultSumw2(1);

  std::string dir = "/afs/cern.ch/work/j/jlawhorn/public/comb_ntuples/";
  double sigscale = 1000;
  double sigscale1 = 10; 
  std::stringstream scale; scale << sigscale;
  std::stringstream scale1; scale1 << sigscale1;

  //Cut definitions
  double luminosity = 3000;
  std::stringstream lumi; lumi << luminosity;

  std::string ttcut = "(tauCat1==1 && tauCat2==1 && bTag1>0 && bTag2>0 && ptTau1>45 && ptTau2>45)";
  std::string mtcut = "((tauCat1==1 && tauCat2==3)||(tauCat2==1 && tauCat1==3))*(bTag1>0 && bTag2>0 && ptTau1>30 && ptTau2>30)";
  std::string etcut = "((tauCat1==1 && tauCat2==2)||(tauCat2==1 && tauCat1==2))*(bTag1>0 && bTag2>0 && ptTau1>30 && ptTau2>30)";
  std::string emcut = "((tauCat1==3 && tauCat2==2)||(tauCat2==3 && tauCat1==2))*(bTag1>0 && bTag2>0 && ptTau1>30 && ptTau2>30)";

  std::string taucut = "(abs(etaTau1)<4.0 && abs(etaTau2)<4.0)";
  std::string bcut = "(ptB1>30 && ptB2>30 && abs(etaB1)<4.0 && abs(etaB2)<4.0)";

  std::string optcut = cut+"*"+taucut+"*"+bcut;
  if (chan_cat==0) { //tt
    optcut += "*"+ttcut;
  }
  else if (chan_cat==1) { //mt
    optcut += "*"+mtcut;
  }
  else if (chan_cat==2) { //et
    optcut += "*"+etcut;
  }
  else if (chan_cat==3) { //em
    optcut += "*"+emcut;
  }

  //signal region
  std::string hhcut = optcut+"*"+scale.str()+"*newEvtWeight*(isBBTT==1)*(1)*"+lumi.str();
  std::string hcut = optcut+"*newEvtWeight*(isBBTT==1)*(1)*"+lumi.str();
  std::string ttbarcut = optcut+"*newEvtWeight*(isBBTT==1)*(1)*"+lumi.str();
  std::string zjetcut = optcut+"*newEvtWeight*(isBBTT==1)*(eventType==4)*"+lumi.str();
  std::string ewkcut = optcut+"*newEvtWeight*(isBBTT==1)*(eventType!=4)*"+lumi.str();

  cout << hhcut << endl;
  
  //--------------------------------------------------------------------------
  
  //Get the trees
  TTree *ttbartree = load(dir+"tt.root"); 
  TTree *htree = load(dir+"H.root"); 
  TTree *ewktree = load(dir+"EWK.root"); 
  TTree *sigtree = load(dir+"HHToTTBB_14TeV.root"); 

  //-------------------------------------------------------------------------

  //Get histograms
  TCanvas *canv0 = MakeCanvas("canv", "histograms", 600, 600);
  canv0->cd();
  std::string vardraw;
  TH1F *sig = new TH1F("HH","",nbins,xmin,xmax);
  vardraw = var+">>"+"HH";
  sigtree->Draw(vardraw.c_str(),hhcut.c_str());
  InitSignal(sig);
  sig->SetLineColor(kBlack);
  TH1F *ttbar = new TH1F("TTbar","",nbins,xmin,xmax);
  vardraw = var+">>"+"TTbar";
  ttbartree->Draw(vardraw.c_str(),ttbarcut.c_str());
  InitHist(ttbar, xtitle.c_str(), ytitle.c_str(), TColor::GetColor(155,152,204), 1001);
  TH1F *hbg = new TH1F("H","",nbins,xmin,xmax);
  vardraw = var+">>"+"H";
  htree->Draw(vardraw.c_str(),hcut.c_str());
  InitHist(hbg, xtitle.c_str(), ytitle.c_str(), TColor::GetColor(141,201,159), 1001);
  TH1F *zjet = new TH1F("Zjets","",nbins,xmin,xmax);
  vardraw = var+">>"+"Zjets";
  ewktree->Draw(vardraw.c_str(),zjetcut.c_str());
  InitHist(zjet, xtitle.c_str(), ytitle.c_str(),  TColor::GetColor(222,90,106), 1001);
  TH1F *ewk = new TH1F("EWK","",nbins,xmin,xmax);
  vardraw = var+">>"+"EWK";
  ewktree->Draw(vardraw.c_str(),ewkcut.c_str());
  InitHist(ewk, xtitle.c_str(), ytitle.c_str(),  TColor::GetColor(248,206,104), 1001);

  cout << sig->GetEntries() << endl;
  cout << ttbar->GetEntries() << endl;
  cout << hbg->GetEntries() << endl;
  cout << zjet->GetEntries() << endl;
  cout << ewk->GetEntries() << endl;
  
  delete canv0;
  //-----------------------------------------------------------------------
  //Draw the histograms
  TCanvas *canv = MakeCanvas("canv", "histograms", 600, 600);
  canv->cd();
  zjet->Add(ewk); hbg->Add(zjet); ttbar->Add(hbg);
  //Error band stat
  TH1F* errorBand = (TH1F*)sig ->Clone("errorBand");
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
  ttbar->SetMaximum(1.2*std::max(maximum(ttbar, 0), maximum(sig, 0)));
  //blind(sig,75,150);
  ttbar->Draw("hist");
  hbg->Draw("histsame");
  zjet->Draw("histsame");
  ewk->Draw("histsame");
  sig->Draw("histsame");
  //errorBand->Draw("e2same");
  canv->RedrawAxis();
  //---------------------------------------------------------------------------
  //Adding a legend
  TLegend* leg = new TLegend(0.53, 0.73, 0.95, 0.90);
  SetLegendStyle(leg);
  leg->AddEntry(sig , TString::Format("%.0f#timeshh#rightarrow#tau#tau bb", sigscale)      , "L");
  leg->AddEntry(ttbar, "t#bar{t}"              , "F" );
  leg->AddEntry(hbg  , "H"  , "F" );
  leg->AddEntry(zjet  , "Z+jets"           , "F" );
  leg->AddEntry(ewk  , "EWK"           , "F" );
  //leg->AddEntry(errorBand,"bkg. uncertainty","F");
  leg->Draw();
  //---------------------------------------------------------------------------
  
  //CMS preliminary 
  char category[50];
  if (chan_cat==0) sprintf(category,"#tau_{h}#tau_{h}");
  if (chan_cat==1) sprintf(category,"#mu#tau_{h}");
  if (chan_cat==2) sprintf(category,"e#tau_{h}");
  if (chan_cat==3) sprintf(category,"e#mu");

  const char* dataset = "CMS Preliminary, HH#rightarrow bb#tau#tau, 3.0 ab^{-1} at 14 TeV";

  //CMSPrelim(dataset, "#tau_{h}#tau_{h}", 0.17, 0.835);
  //CMSPrelim(dataset, "#tau_{h}#tau_{h}", 0.17, 0.835);
  CMSPrelim(dataset, "", 0.16, 0.835);
  TPaveText* chan     = new TPaveText(0.26, 0.77, 0.37, 0.85, "tlbrNDC");
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

  char outfile[50];
  if (chan_cat==0) sprintf(outfile, "%s_tt.png", var.c_str());
  if (chan_cat==1) sprintf(outfile, "%s_mt.png", var.c_str());
  if (chan_cat==2) sprintf(outfile, "%s_et.png", var.c_str());
  if (chan_cat==3) sprintf(outfile, "%s_em.png", var.c_str());

  canv->Print(outfile);

}
