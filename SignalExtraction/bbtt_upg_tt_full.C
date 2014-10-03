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
#include "../Utils/CMS_lumi_v2.h"

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

void bbtt_upg_tt_full(std::string var,int nbins, double xmin, double xmax,std::string xtitle, std::string ytitle, double sigscale=1)
{

  TFile *outDC = new TFile("hh_tt_inputs.root","RECREATE");

  //SetStyle(); gStyle->SetLineStyleString(11,"20 10");
  setTDRStyle();
  TH1::SetDefaultSumw2(1);

  std::string dir = "/data/blue/Bacon/029a/tautaunew/mhhvar/";
  
  std::stringstream scale; scale << sigscale;
  
  //Cut definitions
  double luminosity = 19172;
  std::stringstream lumi; lumi << luminosity;
  std::string objcut = "(pt_1>45 && pt_2>45 && jpt_1>30 && jpt_2>30)*(abs(jeta_1)<2.5 && abs(jeta_2)<2.5 && abs(eta_1)<2.1 && abs(eta_2)<2.1 && byCombinedIsolationDeltaBetaCorrRaw3Hits_2<1.5 &&  byCombinedIsolationDeltaBetaCorrRaw3Hits_1<1.5 && jcsv_1>0.68 && jcsv_2>0.68)";
  //std::string jetcut = objcut+"*(m_sv>100 && m_sv<130 && mjj>90 && mjj<140)";
  std::string jetcut = objcut+"*(1.0)";
  std::string datacut = jetcut+"*(q_1*q_2<0)"; 
  //signal region
  std::string smhhcut = jetcut+"*(weight*decaymodeweight*0.71442e-3/1.972e09)*(q_1*q_2<0)*"+lumi.str();
  std::string tthcut = jetcut+"*(weight*decaymodeweight*0.0789)*(q_1*q_2<0)*"+lumi.str();
  std::string vbfcut = jetcut+"*(weight*decaymodeweight*0.0997)*(q_1*q_2<0)*"+lumi.str();
  std::string gfcut = jetcut+"*(weight*decaymodeweight*1.2179)*(q_1*q_2<0)*"+lumi.str();
  std::string zttmccut = jetcut+"*weight*(decaymodeweight)*(genmatch==3)*(q_1*q_2<0)*"+lumi.str();
  std::string fakecut = jetcut+"*1.06*(q_1*q_2>0)";
  std::string wjetcut = jetcut+"*(q_1*q_2<0)*(weight*decaymodeweight)*"+lumi.str();
  std::string ewkcut = jetcut+"*(q_1*q_2<0)*(weight*decaymodeweight)*"+lumi.str();
  std::string ttbarcut = jetcut+"*1.11*0.96*(q_1*q_2<0)*(weight*decaymodeweight)*"+lumi.str();
  // QCD SS control region
  std::string ttbarcutS = jetcut+"*1.11*0.96*(q_1*q_2>0)*(weight*decaymodeweight)*"+lumi.str();
  //--------------------------------------------------------------------------
  
  //Get the trees
  TTree *datatree = load(dir+"data_select.root"); 
  TTree *hhtree = load(dir+"s12-hh125-bbtt-tt-8tev_select.root");
  TTree *tttree = load(dir+"ttbar-8TeV_select.root");
  TTree *vbfhtree = load(dir+"htt_vtth_sm_125_select.root");
  TTree *gfhtree = load(dir+"htt_gf_sm_125_select.root");
  TTree *assohtree = load(dir+"htt_vtth_sm_125_select.root");
  TTree *vjettree = load(dir+"wjets_select.root");
  TTree *ewktree = load(dir+"ewk-8TeV_select.root");
  TTree *zttmctree = load(dir+"ztt-mad_select.root");
  
  //-------------------------------------------------------------------------
  
  //Get histograms
  TCanvas *canv0 = MakeCanvas("canv", "histograms", 600, 600);
  canv0->cd();
  std::string vardraw;
  TH1F *data = new TH1F("Data","",nbins,xmin,xmax);
  vardraw = var+">>"+"Data";
  datatree->Draw(vardraw.c_str(),datacut.c_str());
  InitHist(data, xtitle.c_str(), ytitle.c_str()); InitData(data);
  TH1F *Ztt = new TH1F("DY","",nbins,xmin,xmax);
  vardraw = var+">>"+"DY";
  zttmctree->Draw(vardraw.c_str(),zttmccut.c_str());
  InitHist(Ztt  , xtitle.c_str(), ytitle.c_str(), TColor::GetColor(248,206,104), 1001);
  TH1F *ttbar = new TH1F("TTbar","",nbins,xmin,xmax);
  vardraw = var+">>"+"TTbar";
  tttree->Draw(vardraw.c_str(),ttbarcut.c_str());
  InitHist(ttbar, xtitle.c_str(), ytitle.c_str(), TColor::GetColor(155,152,204), 1001);
  TH1F *ttbarS = new TH1F("TTbarS","",nbins,xmin,xmax);
  vardraw = var+">>"+"TTbarS";
  tttree->Draw(vardraw.c_str(),ttbarcutS.c_str());
  InitHist(ttbarS, xtitle.c_str(), ytitle.c_str(), TColor::GetColor(155,152,204), 1001);
  TH1F *wjets = new TH1F("Wjets","",nbins,xmin,xmax);
  vardraw = var+">>"+"Wjets";
  vjettree->Draw(vardraw.c_str(),wjetcut.c_str());
  InitHist(wjets, xtitle.c_str(), ytitle.c_str(), TColor::GetColor(222,90,106), 1001);
  TH1F *ewk = new TH1F("Ewk","",nbins,xmin,xmax);
  vardraw = var+">>"+"Ewk";
  ewktree->Draw(vardraw.c_str(),ewkcut.c_str());
  InitHist(ewk, xtitle.c_str(), ytitle.c_str(),  TColor::GetColor(222,90,106), 1001);
  TH1F *qcd = new TH1F("Qcd","",nbins,xmin,xmax);
  vardraw = var+">>"+"Qcd";
  datatree->Draw(vardraw.c_str(),fakecut.c_str());
  InitHist(qcd, xtitle.c_str(), ytitle.c_str(),  kGreen, 1001);
  qcd->Add(ttbarS,-1);
  TH1F *vbfh = new TH1F("VBFH","",nbins,xmin,xmax);
  vardraw = var+">>"+"VBFH";
  vbfhtree->Draw(vardraw.c_str(),vbfcut.c_str());
  InitHist(vbfh, xtitle.c_str(), ytitle.c_str(),  TColor::GetColor(250,202,255), 1001);
  TH1F *ggh = new TH1F("GGH","",nbins,xmin,xmax);
  vardraw = var+">>"+"GGH";
  gfhtree->Draw(vardraw.c_str(),gfcut.c_str());
  InitHist(ggh, xtitle.c_str(), ytitle.c_str(),  TColor::GetColor(250,202,255), 1001);
  TH1F *assoh = new TH1F("AH","",nbins,xmin,xmax);
  vardraw = var+">>"+"AH";
  assohtree->Draw(vardraw.c_str(),tthcut.c_str());
  InitHist(assoh, xtitle.c_str(), ytitle.c_str(),  TColor::GetColor(250,202,255), 1001);
  TH1F *smhh = new TH1F("SMhh","",nbins,xmin,xmax);
  vardraw = var+">>"+"SMhh";
  hhtree->Draw(vardraw.c_str(),smhhcut.c_str());
  InitSignal(smhh);
  smhh->SetLineColor(kBlack); 
  delete canv0;
  //---------------------------------------------------------------------------
  //Print out the yields
  Double_t error=999;
  Double_t sigN=0;
  Double_t sigSig=0;
  Double_t bgdN=0;
  Double_t bgdSig=0;
  //ofstream outfile;
  //outfile.open("yields.txt");
  //outfile << "Yields for the signal region." << std::endl;
  cout << jetcut << endl;
  cout << "SM hh   "  << smhh->IntegralAndError(0,smhh->GetNbinsX(),error) << "+/-";
  sigN+=smhh->IntegralAndError(0,smhh->GetNbinsX(),error);
  sigSig+=error;
  cout << error << endl; error=999;
  cout << "Data   "  << data->IntegralAndError(0,data->GetNbinsX(),error) << "+/-";
  cout << error << endl; error=999;
  cout << "ttbar   "  << ttbar->IntegralAndError(0,ttbar->GetNbinsX(),error) << "+/-";
  cout << error << endl; error=999;
  bgdN+=ttbar->IntegralAndError(0,ttbar->GetNbinsX(),error);bgdSig+=error*error;
  cout << "Ztt     "  << Ztt->IntegralAndError(0,Ztt->GetNbinsX(),error) << "+/-";
  cout << error << endl; error=999;
  bgdN+=Ztt->IntegralAndError(0,Ztt->GetNbinsX(),error);bgdSig+=error*error;
  cout << "ewk     "  << ewk->IntegralAndError(0,ewk->GetNbinsX(),error) << "+/-";
  cout << error << endl; error=999;
  bgdN+=ewk->IntegralAndError(0,ewk->GetNbinsX(),error);bgdSig+=error*error;
  cout << "qcd     "  << qcd->IntegralAndError(0,qcd->GetNbinsX(),error) << "+/-";
  cout << error << endl; error=999;
  bgdN+=qcd->IntegralAndError(0,qcd->GetNbinsX(),error);bgdSig+=error*error;
  cout << "wjets   "  << wjets->IntegralAndError(0,wjets->GetNbinsX(),error) << "+/-";
  cout << error << endl; error=999;
  bgdN+=wjets->IntegralAndError(0,wjets->GetNbinsX(),error);bgdSig+=error*error;
  cout << "ggH     "  << ggh->IntegralAndError(0,ggh->GetNbinsX(),error) << "+/-";
  cout << error << endl; error=999;
  bgdN+=ggh->IntegralAndError(0,ggh->GetNbinsX(),error);bgdSig+=error*error;
  cout << "vbfH    "  << vbfh->IntegralAndError(0,vbfh->GetNbinsX(),error) << "+/-";
  cout << error << endl; error=999;
  bgdN+=vbfh->IntegralAndError(0,vbfh->GetNbinsX(),error);bgdSig+=error*error;
  cout << "vH/ttH  "  << assoh->IntegralAndError(0,assoh->GetNbinsX(),error) << "+/-";
  cout << error << endl; error=999;
  bgdN+=assoh->IntegralAndError(0,assoh->GetNbinsX(),error); bgdSig+=error*error;
  cout << "S = " << sigN << "+/-" << sigSig << endl;
  cout << "B = " << bgdN << "+/-" << TMath::Sqrt(bgdSig) << endl;
  cout << "S/sqrt(B) = " << sigN/TMath::Sqrt(bgdN) << endl;

  //--------------------------------------------------------------------------
  //outfile.close();
  outDC->cd();
  TDirectory* lTD = outDC->mkdir("mutau");
  outDC->cd(lTD->GetPath());
  data->SetName("data_obs");
  data->SetTitle("data_obs");
  data->Write();
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
  outDC->Close();
  //stack some  histtograms together
  vbfh->Add(ggh);
  vbfh->Add(assoh);
  cout << "Single H    "  << vbfh->IntegralAndError(0,vbfh->GetNbinsX(),error) << "+/-";
  cout << error << endl; error=999;
  wjets->Add(ewk); 
  cout << "Total  Electroweak    "  << wjets->IntegralAndError(0,wjets->GetNbinsX(),error) << "+/-";
  cout << error << endl; error=999;
  //-----------------------------------------------------------------------
  smhh->Scale(sigscale);
  //Draw the histograms
  TCanvas *canv = MakeCanvas("canv", "histograms", 800, 600);
  canv->cd();
  ttbar->Add(qcd);
  wjets->Add(ttbar);
  Ztt->Add(wjets);
  vbfh->Add(Ztt);
  //Error band stat
  TH1F* errorBand = (TH1F*)vbfh->Clone("errorBand");
  errorBand  ->SetMarkerSize(0);
  errorBand  ->SetFillColor(13);
  errorBand  ->SetFillStyle(3013);
  errorBand  ->SetLineWidth(1);
 
  data->SetMaximum(2.0*std::max(maximum(data, 0), maximum(data, 0)));
  //blind(data,100,500);
  data->Draw("e");
  vbfh->Draw("histsame");
  Ztt->Draw("histsame");
  wjets->Draw("histsame");
  //ggh->Draw("histsame");
  ttbar->Draw("histsame");
  qcd->Draw("histsame");
  data->Draw("esame");
  errorBand->Draw("e2same");
  smhh->Draw("histsame");
  canv->RedrawAxis();
  //canv->SetLogy(1);
  //---------------------------------------------------------------------------
  //Adding a legend
  TLegend* leg = new TLegend(0.65, 0.65, 0.95, 0.90);
  SetLegendStyle(leg);
  leg->AddEntry(smhh  , TString::Format("%.0f#timeshh#rightarrow#tau#tau bb", sigscale) , "L" );
  leg->AddEntry(data , "Observed"                       , "LP");
  leg->AddEntry(Ztt  , "Z#rightarrow#tau#tau"           , "F" );
  leg->AddEntry(ttbar, "t#bar{t}"                       , "F" );
  leg->AddEntry(wjets  , "Electroweak"                    , "F" );
  leg->AddEntry(vbfh  , "SM H#rightarrow#tau#tau"   , "F" );
  leg->AddEntry(qcd  , "QCD Multi Jet"   , "F" );
  leg->AddEntry(errorBand,"bkg. uncertainty","F");
  leg->Draw();
  //---------------------------------------------------------------------------
   
  CMS_lumi_v2( canv, 8, 11 );
  //CMS preliminary 
  //const char* dataset = "CMS Simulation, 3000 fb^{-1} at 14 TeV";
  const char* category = "";
  //CMSPrelim(dataset, "#tau_{#mu}#tau_{h}", 0.17, 0.835);
  //CMSPrelim(dataset, "", 0.16, 0.835);
  TPaveText* chan     = new TPaveText(0.52, 0.35, 0.91, 0.55, "tlbrNDC");
  chan->SetBorderSize(   0 );
  chan->SetFillStyle(    0 );
  chan->SetTextAlign(   12 );
  chan->SetTextSize ( 0.05 );
  chan->SetTextColor(    1 );
  chan->SetTextFont (   62 );
  chan->AddText(category);
  //chan->Draw();
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
