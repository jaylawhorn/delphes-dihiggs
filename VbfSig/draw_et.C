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


void draw_et(std::string var="mTT",int nbins=10, double xmin=0, double xmax=200,std::string xtitle="m_{#tau#tau}", std::string ytitle="Events")
{
  SetStyle(); gStyle->SetLineStyleString(11,"20 10");
  TH1::SetDefaultSumw2(1);

  std::string dir = "/afs/cern.ch/work/j/jlawhorn/public/vbf/ntuples/";
  double sigscale = 10;
  double sigscale1 = 10; 
  std::stringstream scale; scale << sigscale;
  std::stringstream scale1; scale1 << sigscale1;

  //Cut definitions
  double luminosity = 3000;
  std::stringstream lumi; lumi << luminosity;
  std::string objcut = "(ptTau1>30 && ptTau2>30 && abs(etaTau1) <4.0 && abs(etaTau1)<4.0 )";
  std::string jetcut = objcut+"*(ptJet1>30 && ptJet2>30 && abs(etaJet1) <4.7 && abs(etaJet2) <4.7 )";
  std::string ththcut = jetcut+"*( ( (tauCat1==1 && tauCat2==2) || (tauCat2==1 && tauCat1==2) ) && mJJ>500 && abs(dEta)>3.5)";

  //signal region
  std::string vbfcut = ththcut+"*eventWeight*(eventType==0)*"+lumi.str();
  std::string zttcut = ththcut+"*eventWeight*(eventType==1)*"+lumi.str();
  std::string ttbarcut = ththcut+"*eventWeight*(eventType==3)*"+lumi.str();
  //std::string ewkcut = ththcut+"*eventWeight*(eventType==2)*"+lumi.str();
  std::string othercut = ththcut+"*eventWeight*(eventType==4 || eventType==2)*"+lumi.str();
  
  //--------------------------------------------------------------------------
  
  //Get the trees
  TTree *tree = load(dir+"vbf_comb.root"); 

  //-------------------------------------------------------------------------
  
  //Get histograms
  TCanvas *canv0 = MakeCanvas("canv", "histograms", 600, 600);
  canv0->cd();
  std::string vardraw;
  TH1F *vbf = new TH1F("VBFH","",nbins,xmin,xmax);
  vardraw = var+">>"+"VBFH";
  tree->Draw(vardraw.c_str(),vbfcut.c_str());
  InitSignal(vbf);
  vbf->SetLineColor(kBlack);
  TH1F *ttbar = new TH1F("TTbar","",nbins,xmin,xmax);
  vardraw = var+">>"+"TTbar";
  tree->Draw(vardraw.c_str(),ttbarcut.c_str());
  InitHist(ttbar, xtitle.c_str(), ytitle.c_str(), TColor::GetColor(155,152,204), 1001);
  TH1F *ztt = new TH1F("Ztt","",nbins,xmin,xmax);
  vardraw = var+">>"+"Ztt";
  tree->Draw(vardraw.c_str(),zttcut.c_str());
  InitHist(ztt, xtitle.c_str(), ytitle.c_str(), TColor::GetColor(248,206,104), 1001);
  //TH1F *ewk = new TH1F("Ewk","",nbins,xmin,xmax);
  //vardraw = var+">>"+"Ewk";
  //tree->Draw(vardraw.c_str(),ewkcut.c_str());
  //InitHist(ewk, xtitle.c_str(), ytitle.c_str(),  TColor::GetColor(222,90,106), 1001);
  TH1F *other = new TH1F("Other","",nbins,xmin,xmax);
  vardraw = var+">>"+"Other";
  tree->Draw(vardraw.c_str(),othercut.c_str());
  InitHist(other, xtitle.c_str(), ytitle.c_str(),  TColor::GetColor(222,90,106), 1001);
  
  delete canv0;
  //----------------------------------------------------------------------------
  //Print out the yields
  Double_t error=0.0;
  ofstream outfile;
  outfile.open("yields_et.txt");
  outfile << "Yields for the signal region." << std::endl;
  outfile << "VBF   "  << vbf->IntegralAndError(0,vbf->GetNbinsX(),error) << "+/-" << error << endl;
  outfile << "TTbar   "  << ttbar->IntegralAndError(0,ttbar->GetNbinsX(),error) << "+/-" << error << endl;
  outfile << "Ztt    "  << ztt->IntegralAndError(0,ztt->GetNbinsX(),error) << "+/-" << error << endl;
  //outfile << "ewk    "  << ewk->IntegralAndError(0,ewk->GetNbinsX(),error) << "+/-" << error << endl;
  outfile << "other   "  << other->IntegralAndError(0,other->GetNbinsX(),error) << "+/-" << error << endl;
  outfile << "S/sqrt(B)    "  << vbf->Integral()/(other->Integral()+ztt->Integral()+ttbar->Integral()/*+ewk->Integral()*/) << endl;
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
  outfile.close();
  //-----------------------------------------------------------------------
  //Draw the histograms
  TCanvas *canv = MakeCanvas("canv", "histograms", 600, 600);
  canv->cd();
  //ewk->Add(other); ttbar->Add(ewk); 
  //ztt->Add(ttbar); vbf->Add(ztt);
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
  ztt->SetMaximum(1.2*std::max(maximum(ttbar, 0), maximum(other, 0)));
  //blind(vbf,75,150);
  ztt->Draw("hist");
  ttbar->Draw("histsame");
  //ewk->Draw("histsame");
  other->Draw("histsame");
  vbf->Draw("histsame");
  //errorBand->Draw("e2same");
  canv->RedrawAxis();
  //---------------------------------------------------------------------------
  //Adding a legend
  TLegend* leg = new TLegend(0.53, 0.65, 0.95, 0.90);
  SetLegendStyle(leg);
  leg->AddEntry(vbf , "VBF Signal"             , "F");
  leg->AddEntry(ztt  , "Z#rightarrow#tau#tau"  , "F" );
  leg->AddEntry(ttbar, "t#bar{t}"              , "F" );
  //leg->AddEntry(ewk  , "Electroweak"           , "F" );
  leg->AddEntry(other, "Other"                 , "F" );
  //leg->AddEntry(errorBand,"bkg. uncertainty","F");
  leg->Draw();
  //---------------------------------------------------------------------------
  
  //CMS preliminary 
  const char* dataset = "CMS Preliminary, H#rightarrow#tau#tau, 3.0 ab^{-1} at 14 TeV";
  const char* category = "";
  CMSPrelim(dataset, "e#tau_{h}", 0.17, 0.835);
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
  canv->Print((var+"_et.png").c_str());
  
  /*
    Ratio Data over MC
    
    TCanvas *canv1 = MakeCanvas("canv0", "histograms", 600, 400);
    canv1->SetGridx();
    canv1->SetGridy();
    canv1->cd();
  */
  /*   TH1F* model = (TH1F*)ztt ->Clone("model");
       TH1F* test1 = (TH1F*)data->Clone("test1"); 
       for(int ibin=0; ibin<test1->GetNbinsX(); ++ibin){
       //the small value in case of 0 entries in the model is added to prevent the chis2 test from failing
       model->SetBinContent(ibin+1, model->GetBinContent(ibin+1)>0 ? model->GetBinContent(ibin+1)*model->GetBinWidth(ibin+1) : 0.01);
       //model->SetBinError  (ibin+1, CONVERVATIVE_CHI2 ? 0. : model->GetBinError  (ibin+1)*model->GetBinWidth(ibin+1));
       model->SetBinError  (ibin+1, 0);
     test1->SetBinContent(ibin+1, test1->GetBinContent(ibin+1)*test1->GetBinWidth(ibin+1));
     test1->SetBinError  (ibin+1, test1->GetBinError  (ibin+1)*test1->GetBinWidth(ibin+1));
   }
   double chi2prob = test1->Chi2Test      (model,"PUW");        std::cout << "chi2prob:" << chi2prob << std::endl;
   double chi2ndof = test1->Chi2Test      (model,"CHI2/NDFUW"); std::cout << "chi2ndf :" << chi2ndof << std::endl;
   double ksprob   = test1->KolmogorovTest(model);              std::cout << "ksprob  :" << ksprob   << std::endl;
   double ksprobpe = test1->KolmogorovTest(model,"DX");         std::cout << "ksprobpe:" << ksprobpe << std::endl;  

   std::vector<double> edges;
   TH1F* zero = (TH1F*)ttbar->Clone("zero"); zero->Clear();
   TH1F* rat1 = (TH1F*)data->Clone("rat1"); 
   for(int ibin=0; ibin<rat1->GetNbinsX(); ++ibin){
     rat1->SetBinContent(ibin+1, Ztt->GetBinContent(ibin+1)>0 ? data->GetBinContent(ibin+1)/Ztt->GetBinContent(ibin+1) : 0);
     rat1->SetBinError  (ibin+1, Ztt->GetBinContent(ibin+1)>0 ? data->GetBinError  (ibin+1)/Ztt->GetBinContent(ibin+1) : 0);
     zero->SetBinContent(ibin+1, 0.);
     zero->SetBinError  (ibin+1, Ztt->GetBinContent(ibin+1)>0 ? Ztt ->GetBinError  (ibin+1)/Ztt->GetBinContent(ibin+1) : 0);
   }
   for(int ibin=0; ibin<rat1->GetNbinsX(); ++ibin){
     if(rat1->GetBinContent(ibin+1)>0){
       edges.push_back(TMath::Abs(rat1->GetBinContent(ibin+1)-1.)+TMath::Abs(rat1->GetBinError(ibin+1)));
       // catch cases of 0 bins, which would lead to 0-alpha*0-1
       rat1->SetBinContent(ibin+1, rat1->GetBinContent(ibin+1)-1.);
     }
   }
   float range = 0.1;
   std::sort(edges.begin(), edges.end());
   if (edges[edges.size()-2]>0.1) { range = 0.2; }
   if (edges[edges.size()-2]>0.2) { range = 0.5; }
   if (edges[edges.size()-2]>0.5) { range = 1.0; }
   if (edges[edges.size()-2]>1.0) { range = 1.5; }
   if (edges[edges.size()-2]>1.5) { range = 2.0; }
   rat1->SetLineColor(kBlack);
   rat1->SetFillColor(kGray );
   rat1->SetMaximum(+range);
   rat1->SetMinimum(-range);
   rat1->GetYaxis()->CenterTitle();
   rat1->GetYaxis()->SetTitle("#bf{Data/MC-1}");
   rat1->GetXaxis()->SetTitle("#bf{m_{#tau#tau} [GeV]}");
   rat1->Draw();
   zero->SetFillStyle(  3013);
   zero->SetFillColor(kBlack);
   zero->SetLineColor(kBlack);
   zero->SetMarkerSize(0.1);
   zero->Draw("e2histsame");
   canv1->RedrawAxis();

   TPaveText* stat1 = new TPaveText(0.20, 0.76+0.061, 0.32, 0.76+0.161, "NDC");
   stat1->SetBorderSize(   0 );
   stat1->SetFillStyle(    0 );
   stat1->SetTextAlign(   12 );
   stat1->SetTextSize ( 0.05 );
   stat1->SetTextColor(    1 );
   stat1->SetTextFont (   62 );
   stat1->AddText(TString::Format("#chi^{2}/ndf=%.3f,  P(#chi^{2})=%.3f", chi2ndof, chi2prob));
   //stat1->AddText(TString::Format("#chi^{2}/ndf=%.3f,  P(#chi^{2})=%.3f, P(KS)=%.3f", chi2ndof, chi2prob, ksprob));
   //stat1->Draw();
   canv1->Print((var+"_ratio.png").c_str());*/
}
