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

#include "../Utils/tdrstyle.h"

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


void bbtt_upg_tt(std::string var,int nbins, double xmin, double xmax,std::string xtitle, std::string ytitle, double sigscale=1)
{

  TFile *outDC = new TFile("hh_tt_inputs.root","RECREATE");

  SetStyle(); gStyle->SetLineStyleString(11,"20 10");
  TH1::SetDefaultSumw2(1);
 
  std::string dir = "/afs/cern.ch/work/j/jlawhorn/public/ntuples/";
  
  std::stringstream scale; scale << sigscale;
  
  //Cut definitions
  double luminosity = 3000;
  std::stringstream lumi; lumi << luminosity;
  std::string objcut = "(tauCat1==1 && tauCat2==1 && ptTau1>45 && ptTau2>45 && ptB1>30 && ptB2>30 && (bTag1==2||bTag1==3||bTag1==6||bTag1==7) && (bTag2==1||bTag2==3||bTag2==6||bTag2==7))";
  std::string jetcut = objcut+"*(mTT>50&&mTT<130)*(ptBB1>150)";
  //std::string jetcut = objcut+"*(ptTrk1>0 && ptTrk2>0)*(mTT>50&&mTT<130)*(mBB1>80&&mBB1<140)*(ptBB1>150)*(mHH>300)";
  //signal region
  std::string mccut = jetcut+"*eventWeight*"+lumi.str();
  std::string vbfcut = jetcut+"*eventWeight*49470*0.0632*"+lumi.str();
  std::string sigcut = jetcut+"*eventWeight*"+lumi.str();
  std::string zjetcut = jetcut+"*eventWeight*(eventType!=3&&eventType!=1)*"+lumi.str();
  std::string wjetcut = jetcut+"*eventWeight*(eventType==3&&eventType!=1)*"+lumi.str();
  std::string ewkcut = jetcut+"*eventWeight*(eventType!=1)*"+lumi.str();

  std::string sigcutS = sigcut+"*(ptTrk1>0 && ptTrk2>0)*(mBB1>80&&mBB1<140)";
  std::string mccutS = mccut+"*(ptTrk1>0 && ptTrk2>0)*(mBB1>80&&mBB1<140)";
  std::string vbfcutS = vbfcut+"*(ptTrk1>0 && ptTrk2>0)*(mBB1>80&&mBB1<140)";
  std::string zjetcutS = zjetcut+"*(ptTrk1>0 && ptTrk2>0)*(mBB1>80&&mBB1<140)";
  std::string wjetcutS = wjetcut+"*(ptTrk1>0 && ptTrk2>0)*(mBB1>80&&mBB1<140)";
  std::string ewkcutS = ewkcut+"*(ptTrk1>0 && ptTrk2>0)*(mBB1>80&&mBB1<140)";
  //--------------------------------------------------------------------------
  
  //Get the trees
  TTree *hhtree = load(dir+"HHToTTBB_14TeV.root"); 
  TTree *hhtree_m1 = load(dir+"gFHHTobbtautaulam1m.root");
  TTree *hhtree_m5 = load(dir+"gFHHTobbtautaulam5m.root");
  TTree *hhtree_0 = load(dir+"gFHHTobbtautaulam0.root");
  TTree *hhtree_p5 = load(dir+"gFHHTobbtautaulam5p.root");
  TTree *tttree = load(dir+"tt.root"); 
  TTree *vbfhtree = load(dir+"VBFHToTauTau.root");
  TTree *gfhtree = load(dir+"ggFHToTauTau.root");
  TTree *assohtree = load(dir+"vH_ttH.root");
  TTree *vjettree = load(dir+"Bjets.root");
  TTree *ewktree = load(dir+"ewk.root");
  TTree *qcdtree = load(dir+"qcd.root");
  
  //-------------------------------------------------------------------------
  
  //Get histograms
  TCanvas *canv0 = MakeCanvas("canv", "histograms", 600, 600);
  canv0->cd();
  std::string vardraw;
  Float_t scaleD=1;
  TH1F *Ztt = new TH1F("DY","",nbins,xmin,xmax);
  vardraw = var+">>"+"DY";
  vjettree->Draw(vardraw.c_str(),zjetcut.c_str());
  InitHist(Ztt  , xtitle.c_str(), ytitle.c_str(), TColor::GetColor(248,206,104), 1001);
  TH1F *ZttS = new TH1F("DYS","",nbins,xmin,xmax);
  vardraw = var+">>"+"DYS";
  vjettree->Draw(vardraw.c_str(),zjetcutS.c_str());
  scaleD=ZttS->Integral(0,ZttS->GetNbinsX())/Ztt->Integral(0,Ztt->GetNbinsX());
  Ztt->Scale(scaleD);
  TH1F *ttbar = new TH1F("TTbar","",nbins,xmin,xmax);
  vardraw = var+">>"+"TTbar";
  tttree->Draw(vardraw.c_str(),mccut.c_str());
  InitHist(ttbar, xtitle.c_str(), ytitle.c_str(), TColor::GetColor(155,152,204), 1001);
  TH1F *ttbarS = new TH1F("TTbarS","",nbins,xmin,xmax);
  vardraw = var+">>"+"TTbarS";
  tttree->Draw(vardraw.c_str(),mccutS.c_str());
  scaleD=ttbarS->Integral(0,ttbarS->GetNbinsX())/ttbar->Integral(0,ttbar->GetNbinsX());
  ttbar->Scale(scaleD);
  TH1F *wjets = new TH1F("Wjets","",nbins,xmin,xmax);
  vardraw = var+">>"+"Wjets";
  vjettree->Draw(vardraw.c_str(),wjetcut.c_str());
  InitHist(wjets, xtitle.c_str(), ytitle.c_str(), TColor::GetColor(222,90,106), 1001);
  TH1F *wjetsS = new TH1F("WjetsS","",nbins,xmin,xmax);
  vardraw = var+">>"+"WjetsS";
  vjettree->Draw(vardraw.c_str(),wjetcutS.c_str());
  scaleD=wjetsS->Integral(0,wjetsS->GetNbinsX())/wjets->Integral(0,wjets->GetNbinsX());
  wjets->Scale(scaleD);
  TH1F *ewk = new TH1F("Ewk","",nbins,xmin,xmax);
  vardraw = var+">>"+"Ewk";
  ewktree->Draw(vardraw.c_str(),ewkcut.c_str());
  InitHist(ewk, xtitle.c_str(), ytitle.c_str(),  TColor::GetColor(222,90,106), 1001);
  TH1F *ewkS = new TH1F("EwkS","",nbins,xmin,xmax);
  vardraw = var+">>"+"EwkS";
  ewktree->Draw(vardraw.c_str(),ewkcutS.c_str());
  scaleD=ewkS->Integral(0,ewkS->GetNbinsX())/ewk->Integral(0,ewk->GetNbinsX());
  ewk->Scale(scaleD);
  TH1F *qcd = new TH1F("Qcd","",nbins,xmin,xmax);
  vardraw = var+">>"+"Qcd";
  qcdtree->Draw(vardraw.c_str(),mccut.c_str());
  InitHist(qcd, xtitle.c_str(), ytitle.c_str(),  TColor::GetColor(222,90,106), 1001);
  TH1F *qcdS = new TH1F("QcdS","",nbins,xmin,xmax);
  vardraw = var+">>"+"QcdS";
  qcdtree->Draw(vardraw.c_str(),mccutS.c_str());
  scaleD=qcdS->Integral(0,qcdS->GetNbinsX())/qcd->Integral(0,qcd->GetNbinsX());
  qcd->Scale(scaleD);
  TH1F *vbfh = new TH1F("VBFH","",nbins,xmin,xmax);
  vardraw = var+">>"+"VBFH";
  vbfhtree->Draw(vardraw.c_str(),vbfcut.c_str());
  InitHist(vbfh, xtitle.c_str(), ytitle.c_str(),  TColor::GetColor(250,202,255), 1001);
  TH1F *vbfhS = new TH1F("VbfhS","",nbins,xmin,xmax);
  vardraw = var+">>"+"VbfhS";
  vbfhtree->Draw(vardraw.c_str(),vbfcutS.c_str());
  scaleD=vbfhS->Integral(0,vbfhS->GetNbinsX())/vbfh->Integral(0,vbfh->GetNbinsX());
  vbfh->Scale(scaleD);
  TH1F *ggh = new TH1F("GGH","",nbins,xmin,xmax);
  vardraw = var+">>"+"GGH";
  gfhtree->Draw(vardraw.c_str(),mccut.c_str());
  InitHist(ggh, xtitle.c_str(), ytitle.c_str(),  TColor::GetColor(250,202,255), 1001);
  TH1F *gghS = new TH1F("GghS","",nbins,xmin,xmax);
  vardraw = var+">>"+"GghS";
  gfhtree->Draw(vardraw.c_str(),mccutS.c_str());
  scaleD=gghS->Integral(0,gghS->GetNbinsX())/ggh->Integral(0,ggh->GetNbinsX());
  ggh->Scale(scaleD);
  TH1F *assoh = new TH1F("AH","",nbins,xmin,xmax);
  vardraw = var+">>"+"AH";
  assohtree->Draw(vardraw.c_str(),mccut.c_str());
  InitHist(assoh, xtitle.c_str(), ytitle.c_str(),  TColor::GetColor(250,202,255), 1001);
  TH1F *assohS = new TH1F("AssohS","",nbins,xmin,xmax);
  vardraw = var+">>"+"AssohS";
  assohtree->Draw(vardraw.c_str(),mccutS.c_str());
  scaleD=assohS->Integral(0,assohS->GetNbinsX())/assoh->Integral(0,assoh->GetNbinsX());
  assoh->Scale(scaleD);
  TH1F *smhh = new TH1F("SMhh","",nbins,xmin,xmax);
  vardraw = var+">>"+"SMhh";
  hhtree->Draw(vardraw.c_str(),sigcut.c_str());
  InitSignal(smhh);
  smhh->SetLineColor(kBlack);
  TH1F *smhhS = new TH1F("SMhhS","",nbins,xmin,xmax);
  vardraw = var+">>"+"SMhhS";
  hhtree->Draw(vardraw.c_str(),sigcutS.c_str());
  scaleD=smhhS->Integral(0,smhhS->GetNbinsX())/smhh->Integral(0,smhh->GetNbinsX());
  smhh->Scale(scaleD);
  TH1F *hh_0 = new TH1F("hh_0","",nbins,xmin,xmax);
  vardraw = var+">>"+"hh_0";
  hhtree_0->Draw(vardraw.c_str(),sigcut.c_str());
  InitSignal(hh_0);
  hh_0->SetLineColor(kBlack);
  TH1F *hh_m1 = new TH1F("hh_m1","",nbins,xmin,xmax);
  vardraw = var+">>"+"hh_m1";
  hhtree_m1->Draw(vardraw.c_str(),sigcut.c_str());
  InitSignal(hh_m1);
  hh_m1->SetLineColor(kBlack);
  TH1F *hh_m5 = new TH1F("hh_m5","",nbins,xmin,xmax);
  vardraw = var+">>"+"hh_m5";
  hhtree_m5->Draw(vardraw.c_str(),sigcut.c_str());
  InitSignal(hh_m5);
  hh_m5->SetLineColor(kBlack);
  TH1F *hh_p5 = new TH1F("hh_p5","",nbins,xmin,xmax);
  vardraw = var+">>"+"hh_p5";
  hhtree_p5->Draw(vardraw.c_str(),sigcut.c_str());
  InitSignal(hh_p5);
  hh_p5->SetLineColor(kBlack);
  delete canv0;
  //---------------------------------------------------------------------------
  //Print out the yields
  Double_t error=999;
  Double_t sigN=0;
  Double_t sigSig=0;
  Double_t bgdN=0;
  Double_t bgdSig=0;
  Double_t bgdN1=0;
  Double_t bgdSig1=0;
  //ofstream outfile;
  //outfile.open("yields.txt");
  //outfile << "Yields for the signal region." << std::endl;
  cout << jetcut << endl;
  cout << "SM hh   "  << smhh->IntegralAndError(0,smhh->GetNbinsX(),error) << "+/-";
  cout << error << endl; error=999;
  cout << "SM hh S "  << smhhS->IntegralAndError(0,smhhS->GetNbinsX(),error) << "+/-";
  cout << error << endl; error=999;
  sigN+=smhhS->IntegralAndError(0,smhhS->GetNbinsX(),error);
  sigSig+=error;
  cout << " 0 x lam "  << hh_0->IntegralAndError(0,hh_0->GetNbinsX(),error) << "+/-";
  cout << error << endl; error=999;
  cout << "-1 x lam "  << hh_m1->IntegralAndError(0,hh_m1->GetNbinsX(),error) << "+/-";
  cout << error << endl; error=999;
  cout << "-5 x lam "  << hh_m5->IntegralAndError(0,hh_m5->GetNbinsX(),error) << "+/-";
  cout << error << endl; error=999;
  cout << " 5 x lam "  << hh_p5->IntegralAndError(0,hh_p5->GetNbinsX(),error) << "+/-";
  cout << error << endl; error=999;
  cout << "ttbar   "  << ttbar->IntegralAndError(0,ttbar->GetNbinsX(),error) << "+/-";
  cout << error << endl; error=999;
  cout << "ttbar S "  << ttbarS->IntegralAndError(0,ttbarS->GetNbinsX(),error) << "+/-";
  cout << error << endl; error=999;
  bgdN+=ttbar->IntegralAndError(0,ttbar->GetNbinsX(),error); bgdSig+=error*error;
  bgdN1+=ttbarS->IntegralAndError(0,ttbarS->GetNbinsX(),error); bgdSig1+=error*error;
  cout << "Ztt     "  << Ztt->IntegralAndError(0,Ztt->GetNbinsX(),error) << "+/-";
  cout << error << endl; error=999;
  cout << "Ztt S   "  << ZttS->IntegralAndError(0,ZttS->GetNbinsX(),error) << "+/-";
  cout << error << endl; error=999;
  bgdN+=Ztt->IntegralAndError(0,Ztt->GetNbinsX(),error); bgdSig+=error*error;
  bgdN1+=ZttS->IntegralAndError(0,ZttS->GetNbinsX(),error); bgdSig1+=error*error;
  cout << "ewk     "  << ewk->IntegralAndError(0,ewk->GetNbinsX(),error) << "+/-";
  cout << error << endl; error=999;
  cout << "ewk S   "  << ewkS->IntegralAndError(0,ewkS->GetNbinsX(),error) << "+/-";
  cout << error << endl; error=999;
  bgdN+=ewk->IntegralAndError(0,ewk->GetNbinsX(),error); bgdSig+=error*error;
  bgdN1+=ewkS->IntegralAndError(0,ewkS->GetNbinsX(),error); bgdSig1+=error*error;
  cout << "qcd     "  << qcd->IntegralAndError(0,qcd->GetNbinsX(),error) << "+/-";
  cout << error << endl; error=999;
  cout << "qcd S   "  << qcdS->IntegralAndError(0,qcdS->GetNbinsX(),error) << "+/-";
  cout << error << endl; error=999;
  //bgdN+=qcd->IntegralAndError(0,qcd->GetNbinsX(),error);
  cout << "wjets   "  << wjets->IntegralAndError(0,wjets->GetNbinsX(),error) << "+/-";
  cout << error << endl; error=999;
  cout << "wjets S "  << wjetsS->IntegralAndError(0,wjetsS->GetNbinsX(),error) << "+/-";
  cout << error << endl; error=999;
  bgdN+=wjets->IntegralAndError(0,wjets->GetNbinsX(),error); bgdSig+=error*error;
  bgdN1+=wjetsS->IntegralAndError(0,wjetsS->GetNbinsX(),error); bgdSig1+=error*error;
  cout << "ggH     "  << ggh->IntegralAndError(0,ggh->GetNbinsX(),error) << "+/-";
  cout << error << endl; error=999;
  cout << "ggH S   "  << gghS->IntegralAndError(0,gghS->GetNbinsX(),error) << "+/-";
  cout << error << endl; error=999;
  bgdN+=ggh->IntegralAndError(0,ggh->GetNbinsX(),error); bgdSig+=error*error;
  bgdN1+=gghS->IntegralAndError(0,gghS->GetNbinsX(),error); bgdSig1+=error*error;
  cout << "vbfH    "  << vbfh->IntegralAndError(0,vbfh->GetNbinsX(),error) << "+/-";
  cout << error << endl; error=999;
  cout << "vbfH S  "  << vbfhS->IntegralAndError(0,vbfhS->GetNbinsX(),error) << "+/-";
  cout << error << endl; error=999;
  bgdN+=vbfh->IntegralAndError(0,vbfh->GetNbinsX(),error); bgdSig+=error*error;
  bgdN1+=vbfhS->IntegralAndError(0,vbfhS->GetNbinsX(),error); bgdSig1+=error*error;
  cout << "vH/ttH  "  << assoh->IntegralAndError(0,assoh->GetNbinsX(),error) << "+/-";
  cout << error << endl; error=999;
  cout << "v/ttH S "  << assohS->IntegralAndError(0,assohS->GetNbinsX(),error) << "+/-";
  cout << error << endl; error=999;
  bgdN+=assoh->IntegralAndError(0,assoh->GetNbinsX(),error); bgdSig+=error*error;
  bgdN1+=assohS->IntegralAndError(0,assohS->GetNbinsX(),error); bgdSig1+=error*error;
  cout << "S = " << sigN << "+/-" << sigSig << endl;
  cout << "B = " << bgdN << "+/-" << TMath::Sqrt(bgdSig) << endl;
  cout << "B = " << bgdN1 << "+/-" << TMath::Sqrt(bgdSig1) << endl;
  cout << "S/sqrt(B) = " << sigN/TMath::Sqrt(bgdN) << endl;

  //--------------------------------------------------------------------------
  //outfile.close();
  outDC->cd();
  TDirectory* lTD = outDC->mkdir("tautau");
  outDC->cd(lTD->GetPath());
  ttbar->SetName("data_obs");
  ttbar->SetTitle("data_obs");
  ttbar->Write();
  Ztt->SetName("ZTT");
  Ztt->SetTitle("ZTT");
  Ztt->Write();
  ttbar->SetName("TT");
  ttbar->SetTitle("TT");
  ttbar->Write();
  wjets->SetName("W");
  wjets->SetTitle("W");
  wjets->Write();
  ewk->SetName("VV");
  ewk->SetTitle("VV");
  ewk->Write();
  ewk->SetName("VV");
  ewk->SetTitle("VV");
  ewk->Write();
  qcd->SetName("QCD");
  qcd->SetTitle("QCD");
  qcd->Write();
  vbfh->SetName("qqH");
  vbfh->SetTitle("qqH");
  vbfh->Write();
  ggh->SetName("ggH");
  ggh->SetTitle("ggH");
  ggh->Write();
  assoh->SetName("assoH");
  assoh->SetTitle("assoH");
  assoh->Write();
  smhh->SetName("ggHH");
  smhh->SetTitle("ggHH");
  smhh->Write();
  hh_0->SetName("lam0");
  hh_0->SetTitle("lam0");
  hh_0->Write();
  hh_m5->SetName("lamm5");
  hh_m5->SetTitle("lamm5");
  hh_m5->Write();
  hh_m1->SetName("lamm1");
  hh_m1->SetTitle("lamm1");
  hh_m1->Write();
  hh_p5->SetName("lamp5");
  hh_p5->SetTitle("lamp5");
  hh_p5->Write();
  outDC->Close();
  //stack some  histtograms together
  vbfh->Add(ggh); 
  vbfh->Add(assoh);
  wjets->Add(ewk); 
  //-----------------------------------------------------------------------
  smhh->Scale(sigscale);
  //Draw the histograms
  TCanvas *canv = MakeCanvas("canv", "histograms", 600, 600);
  canv->cd();
  wjets->Add(ttbar);
  Ztt->Add(wjets);
  vbfh->Add(Ztt);
  //Error band stat
  TH1F* errorBand = (TH1F*)vbfh->Clone("errorBand");
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
  vbfh->SetMaximum(1.1*std::max(maximum(vbfh, 0), maximum(smhh, 0)));
  //blind(data,75,150);
  //data->Draw("e");
  vbfh->Draw("hist");
  //qcd->Draw("histsame");
  Ztt->Draw("histsame");
  wjets->Draw("histsame");
  //ggh->Draw("histsame");
  ttbar->Draw("histsame");
  //data->Draw("esame");
  errorBand->Draw("e2same");
  smhh->Draw("histsame");
  canv->RedrawAxis();
  //canv->SetLogy(1);
  //---------------------------------------------------------------------------
  //Adding a legend
  TLegend* leg = new TLegend(0.65, 0.65, 0.95, 0.90);
  SetLegendStyle(leg);
  leg->AddEntry(smhh  , TString::Format("%.0f#timeshh#rightarrow#tau#tau bb", sigscale) , "L" );
  //leg->AddEntry(data , "Observed"                       , "LP");
  leg->AddEntry(Ztt  , "Z#rightarrow#tau#tau"           , "F" );
  leg->AddEntry(ttbar, "t#bar{t}"                       , "F" );
  leg->AddEntry(wjets  , "Electroweak"                    , "F" );
  leg->AddEntry(vbfh  , "SM H#rightarrow#tau#tau"   , "F" );
  leg->AddEntry(errorBand,"bkg. uncertainty","F");
  leg->Draw();
  //---------------------------------------------------------------------------
   
  //CMS preliminary 
  const char* dataset = "CMS Simulation, 3000 fb^{-1} at 14 TeV";
  const char* category = "";
  CMSPrelim(dataset, "#tau_{h}#tau_{h}", 0.17, 0.835);
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
  canv->Print((var+"_tt.png").c_str());
  
  /*
    Ratio Data over MC
  */
  /*
  TCanvas *canv1 = MakeCanvas("canv0", "histograms", 600, 400);
  canv1->SetGridx();
  canv1->SetGridy();
  canv1->cd();

  TH1F* model = (TH1F*)Ztt ->Clone("model");
  TH1F* test1 = (TH1F*)vbfh->Clone("test1"); 
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
  TH1F* rat1 = (TH1F*)vbfh->Clone("rat1"); 
  for(int ibin=0; ibin<rat1->GetNbinsX(); ++ibin){
    rat1->SetBinContent(ibin+1, Ztt->GetBinContent(ibin+1)>0 ? vbfh->GetBinContent(ibin+1)/Ztt->GetBinContent(ibin+1) : 0);
    rat1->SetBinError  (ibin+1, Ztt->GetBinContent(ibin+1)>0 ? vbfh->GetBinError  (ibin+1)/Ztt->GetBinContent(ibin+1) : 0);
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
  canv1->Print((var+"_ratio.png").c_str());
  */
}
