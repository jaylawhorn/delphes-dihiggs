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

void bbtt_eff_tt(std::string var="mt2pileup",int nbins=10, double xmin=0, double xmax=600,std::string xtitle="", std::string ytitle="", double sigscale=1,int hist=1)
{
  double massLEdges[14]    = {-0.5,-0.45,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.0,0.05,0.1,0.2};
  TFile *outDC = new TFile("hh_tt_inputs.root","RECREATE");

  //SetStyle(); gStyle->SetLineStyleString(11,"20 10");
  //setTDRStyle();
  TH1::SetDefaultSumw2(1);

  //std::string dir = "/data/blue/Bacon/029a/Upgrade/sep-2-bbtt/";
  std::string dir = "/afs/cern.ch/work/j/jlawhorn/para-eff/";
  
  std::stringstream scale; scale << sigscale;

  // Tau eff 65%, fake {(abs(eta)<1.8)*0.006+(abs(eta)>1.8)*0.015}
  // b.... complicated
  double eff_points[10] = { 0.55, 0.57, 0.59, 0.61, 0.63, 0.65, 0.67, 0.69, 0.71, 0.73 };
  
  //Cut definitions
  double luminosity = 3000;
  std::stringstream lumi; lumi << luminosity;

  std::string objcut = "(tauCat1==1 && tauCat2==1 && (ptTau1>45 && ptTau2>45) && ptB1>30 && ptB2>30 && (bTag1==2||bTag1==3||bTag1==6||bTag1==7) && (bTag2==2||bTag2==3||bTag2==6||bTag2==7))*(abs(etaTau1)<2.1 && abs(etaTau2)<2.1 && abs(etaB1)<2.5 && abs(etaB2)<2.5 && ptTrk1>0 && ptTrk2>0)";
  std::string jetcut = objcut+"*((m_svpileup>90 && m_svpileup<120 && (mBB1>90&&mBB1<140)))";

  std::string real2cut = "*(rTau1==1&&rTau2==1)";
  std::string real1cut = "*((rTau1==1&&rTau2==0)||(rTau1==0&&rTau2==1))";
  std::string fakecut = "*(rTau1==0&&rTau2==0)";
  //signal region
  std::string sigreal2cut = jetcut+real2cut+"*eventWeight*"+lumi.str();
  std::string sigreal1cut = jetcut+real1cut+"*eventWeight*"+lumi.str();
  std::string sigfakecut = jetcut+fakecut+"*eventWeight*"+lumi.str();

  std::string mcreal2cut = jetcut+real2cut+"*eventWeight*"+lumi.str();
  std::string mcreal1cut = jetcut+real1cut+"*eventWeight*"+lumi.str();
  std::string mcfakecut = jetcut+fakecut+"*eventWeight*"+lumi.str();

  std::string zjetreal2cut = jetcut+real2cut+"*eventWeight*(eventType!=3&&eventType!=1)*"+lumi.str();
  std::string zjetreal1cut = jetcut+real1cut+"*eventWeight*(eventType!=3&&eventType!=1)*"+lumi.str();
  std::string zjetfakecut = jetcut+fakecut+"*eventWeight*(eventType!=3&&eventType!=1)*"+lumi.str();

  std::string wjetreal2cut = jetcut+real2cut+"*eventWeight*(eventType==3&&eventType!=1)*"+lumi.str();
  std::string wjetreal1cut = jetcut+real1cut+"*eventWeight*(eventType==3&&eventType!=1)*"+lumi.str();
  std::string wjetfakecut = jetcut+fakecut+"*eventWeight*(eventType==3&&eventType!=1)*"+lumi.str();

  std::string ewkreal2cut = jetcut+real2cut+"*eventWeight*(eventType!=1)*"+lumi.str();
  std::string ewkreal1cut = jetcut+real1cut+"*eventWeight*(eventType!=1)*"+lumi.str();
  std::string ewkfakecut = jetcut+fakecut+"*eventWeight*(eventType!=1)*"+lumi.str();

  std::string vbfreal2cut = jetcut+real2cut+"*eventWeight*49470*0.0632*"+lumi.str();  
  std::string vbfreal1cut = jetcut+real1cut+"*eventWeight*49470*0.0632*"+lumi.str();  
  std::string vbffakecut = jetcut+fakecut+"*eventWeight*49470*0.0632*"+lumi.str();  

  //--------------------------------------------------------------------------
  
  //Get the trees
  //TTree *hhtree = load(dir+"HHToTTBB_14TeV-small.root");
  TTree *hhtree = load(dir+"gFHHTobbtautau-small.root");
  //TTree *hhtree_m1 = load(dir+"gFHHTobbtautaulam1m-small.root");
  //TTree *hhtree_m5 = load(dir+"gFHHTobbtautaulam5m-small.root");
  //TTree *hhtree_0 = load(dir+"gFHHTobbtautaulam0-small.root");
  //TTree *hhtree_p5 = load(dir+"gFHHTobbtautaulam5p-small.root");
  TTree *hhtree_vbf = load(dir+"vbfHHTobbtautau-small.root");
  TTree *tttree = load(dir+"tt-small.root");
  TTree *vbfhtree = load(dir+"VBFHToTauTau-small.root");
  TTree *gfhtree = load(dir+"ggFHToTauTau-small.root");
  TTree *assohtree = load(dir+"vH_ttH-small.root");
  TTree *vjettree = load(dir+"Bjets-small.root");
  TTree *ewktree = load(dir+"ewk-small.root");
  
  //-------------------------------------------------------------------------
  
  //Get histograms
  TCanvas *canv0 = MakeCanvas("canv", "histograms", 600, 600);
  canv0->cd();

  std::string vardraw;
  std::string real2name;
  std::string real1name;
  std::string fakename;

  for (Int_t i=0; i<10; i++) {
    stringstream j; j<< i;

    double scaleR=eff_points[i]/0.65;
    stringstream scaleRstr; scaleRstr << scaleR;
    //cout << scaleR << endl;
    
    TH1F *smhh, *smhh_real2, *smhh_real1, *smhh_fake;
    if(hist) {
      smhh= new TH1F("SMhh","",nbins,xmin,xmax);
      real2name="SMhh_real_"+j.str();
      smhh_real2= new TH1F(real2name.c_str(),"",nbins,xmin,xmax);
      real1name="SMhh_half_"+j.str();
      smhh_real1= new TH1F(real1name.c_str(),"",nbins,xmin,xmax);
      fakename="SMhh_fake_"+j.str();
      smhh_fake= new TH1F(fakename.c_str(),"",nbins,xmin,xmax);
    }
    else {
      smhh= new TH1F("SMhh","",13,massLEdges);
      real2name="SMhh_real_"+j.str();
      smhh_real2= new TH1F(real2name.c_str(),"",13,massLEdges);
      real1name="SMhh_half_"+j.str();
      smhh_real1= new TH1F(real1name.c_str(),"",13,massLEdges);
      fakename="SMhh_fake_"+j.str();
      smhh_fake= new TH1F(fakename.c_str(),"",13,massLEdges);
    }
    vardraw = var+">>"+fakename;
    hhtree->Draw(vardraw.c_str(),sigfakecut.c_str());
    vardraw = var+">>"+real1name;
    hhtree->Draw(vardraw.c_str(),sigreal1cut.c_str());
    vardraw = var+">>"+real2name;
    hhtree->Draw(vardraw.c_str(),sigreal2cut.c_str());

    TH1F *Ztt, *Ztt_real2, *Ztt_real1, *Ztt_fake;
    if(hist) {
      Ztt = new TH1F("DY","",nbins,xmin,xmax);
      real2name="DY_real_"+j.str();
      Ztt_real2 = new TH1F(real2name.c_str(),"",nbins,xmin,xmax);
      real1name="DY_half_"+j.str();
      Ztt_real1 = new TH1F(real1name.c_str(),"",nbins,xmin,xmax);
      fakename="DY_fake_"+j.str();
      Ztt_fake = new TH1F(fakename.c_str(),"",nbins,xmin,xmax);
    }
    else {
      Ztt = new TH1F("DY","",13,massLEdges);
      real2name="DY_real_"+j.str();
      Ztt_real2 = new TH1F(real2name.c_str(),"",13,massLEdges);
      real1name="DY_half_"+j.str();
      Ztt_real1 = new TH1F(real1name.c_str(),"",13,massLEdges);
      fakename="DY_fake_"+j.str();
      Ztt_fake = new TH1F(fakename.c_str(),"",13,massLEdges);
    }
    vardraw = var+">>"+fakename;
    vjettree->Draw(vardraw.c_str(),zjetfakecut.c_str());
    vardraw = var+">>"+real1name;
    vjettree->Draw(vardraw.c_str(),zjetreal1cut.c_str());
    vardraw = var+">>"+real2name;
    vjettree->Draw(vardraw.c_str(),zjetreal2cut.c_str());

    TH1F *Wjets, *Wjets_real2, *Wjets_real1, *Wjets_fake;
    if(hist) {
      Wjets = new TH1F("Wjets","",nbins,xmin,xmax);
      real2name="Wjets_real_"+j.str();
      Wjets_real2 = new TH1F(real2name.c_str(),"",nbins,xmin,xmax);
      real1name="Wjets_half_"+j.str();
      Wjets_real1 = new TH1F(real1name.c_str(),"",nbins,xmin,xmax);
      fakename="Wjets_fake_"+j.str();
      Wjets_fake = new TH1F(fakename.c_str(),"",nbins,xmin,xmax);
    }
    else {
      Wjets = new TH1F("Wjets","",13,massLEdges);
      real2name="Wjets_real_"+j.str();
      Wjets_real2 = new TH1F(real2name.c_str(),"",13,massLEdges);
      real1name="Wjets_half_"+j.str();
      Wjets_real1 = new TH1F(real1name.c_str(),"",13,massLEdges);
      fakename="Wjets_fake_"+j.str();
      Wjets_fake = new TH1F(fakename.c_str(),"",13,massLEdges);
    }
    vardraw = var+">>"+fakename;
    vjettree->Draw(vardraw.c_str(),wjetfakecut.c_str());
    vardraw = var+">>"+real2name;
    vjettree->Draw(vardraw.c_str(),wjetreal2cut.c_str());
    vardraw = var+">>"+real1name;
    vjettree->Draw(vardraw.c_str(),wjetreal1cut.c_str());

    TH1F *ewk, *ewk_real2, *ewk_real1, *ewk_fake;
    if(hist) {
      ewk = new TH1F("Ewk","",nbins,xmin,xmax);
      real2name="Ewk_real_"+j.str();
      ewk_real2 = new TH1F(real2name.c_str(),"",nbins,xmin,xmax);
      real1name="Ewk_half_"+j.str();
      ewk_real1 = new TH1F(real1name.c_str(),"",nbins,xmin,xmax);
      fakename="Ewk_fake_"+j.str();
      ewk_fake = new TH1F(fakename.c_str(),"",nbins,xmin,xmax);
    }
    else {
      ewk = new TH1F("Ewk","",13,massLEdges);
      real2name="Ewk_real_"+j.str();
      ewk_real2 = new TH1F(real2name.c_str(),"",13,massLEdges);
      real1name="Ewk_half_"+j.str();
      ewk_real1 = new TH1F(real1name.c_str(),"",13,massLEdges);
      fakename="Ewk_fake_"+j.str();
      ewk_fake = new TH1F(fakename.c_str(),"",13,massLEdges);
    }
    vardraw = var+">>"+fakename;
    ewktree->Draw(vardraw.c_str(),ewkfakecut.c_str());
    vardraw = var+">>"+real2name;
    ewktree->Draw(vardraw.c_str(),ewkreal2cut.c_str());
    vardraw = var+">>"+real1name;
    ewktree->Draw(vardraw.c_str(),ewkreal1cut.c_str());

    TH1F *vbfh, *vbfh_real2, *vbfh_real1, *vbfh_fake;
    if(hist) {
      vbfh = new TH1F("VBFH","",nbins,xmin,xmax);
      real2name="VBFH_real_"+j.str();
      vbfh_real2 = new TH1F(real2name.c_str(),"",nbins,xmin,xmax);
      real1name="VBFH_half_"+j.str();
      vbfh_real1 = new TH1F(real1name.c_str(),"",nbins,xmin,xmax);
      fakename="VBFH_fake_"+j.str();
      vbfh_fake = new TH1F(fakename.c_str(),"",nbins,xmin,xmax);
    }
    else {
      vbfh = new TH1F("VBFH","",13,massLEdges);
      real2name="VBFH_real_"+j.str();
      vbfh_real2 = new TH1F(real2name.c_str(),"",13,massLEdges);
      real1name="VBFH_half_"+j.str();
      vbfh_real1 = new TH1F(real1name.c_str(),"",13,massLEdges);
      fakename="VBFH_fake_"+j.str();
      vbfh_fake = new TH1F(fakename.c_str(),"",13,massLEdges);
    }
    vardraw = var+">>"+fakename;
    vbfhtree->Draw(vardraw.c_str(),vbffakecut.c_str());
    vardraw = var+">>"+real2name;
    vbfhtree->Draw(vardraw.c_str(),vbfreal2cut.c_str());
    vardraw = var+">>"+real1name;
    vbfhtree->Draw(vardraw.c_str(),vbfreal1cut.c_str());

    TH1F *ggh, *ggh_real2, *ggh_real1, *ggh_fake;
    if(hist) {
      ggh = new TH1F("GGH","",nbins,xmin,xmax);
      real2name="GGH_real_"+j.str();
      ggh_real2 = new TH1F(real2name.c_str(),"",nbins,xmin,xmax);
      real1name="GGH_half_"+j.str();
      ggh_real1 = new TH1F(real1name.c_str(),"",nbins,xmin,xmax);
      fakename="GGH_fake_"+j.str();
      ggh_fake = new TH1F(fakename.c_str(),"",nbins,xmin,xmax);
    }
    else {
      ggh = new TH1F("GGH","",13,massLEdges);
      real2name="GGH_real_"+j.str();
      ggh_real2 = new TH1F(real2name.c_str(),"",13,massLEdges);
      real1name="GGH_half_"+j.str();
      ggh_real1 = new TH1F(real1name.c_str(),"",13,massLEdges);
      fakename="GGH_fake_"+j.str();
      ggh_fake = new TH1F(fakename.c_str(),"",13,massLEdges);
    }
    vardraw = var+">>"+fakename;
    gfhtree->Draw(vardraw.c_str(),mcfakecut.c_str());
    vardraw = var+">>"+real2name;
    gfhtree->Draw(vardraw.c_str(),mcreal2cut.c_str());
    vardraw = var+">>"+real1name;
    gfhtree->Draw(vardraw.c_str(),mcreal1cut.c_str());

    TH1F *assoh, *assoh_real2, *assoh_real1, *assoh_fake;
    if(hist) {
      assoh = new TH1F("AH","",nbins,xmin,xmax);
      real2name="AH_real_"+j.str();
      assoh_real2 = new TH1F(real2name.c_str(),"",nbins,xmin,xmax);
      real1name="AH_half_"+j.str();
      assoh_real1 = new TH1F(real1name.c_str(),"",nbins,xmin,xmax);
      fakename="AH_fake_"+j.str();
      assoh_fake = new TH1F(fakename.c_str(),"",nbins,xmin,xmax);
    }
    else {
      assoh = new TH1F("AH","",13,massLEdges);
      real2name="AH_real_"+j.str();
      assoh_real2 = new TH1F(real2name.c_str(),"",13,massLEdges);
      real1name="AH_half_"+j.str();
      assoh_real1 = new TH1F(real1name.c_str(),"",13,massLEdges);
      fakename="AH_fake_"+j.str();
      assoh_fake = new TH1F(fakename.c_str(),"",13,massLEdges);
    }
    vardraw = var+">>"+fakename;
    assohtree->Draw(vardraw.c_str(),mcfakecut.c_str());
    vardraw = var+">>"+real2name;
    assohtree->Draw(vardraw.c_str(),mcreal2cut.c_str());
    vardraw = var+">>"+real1name;
    assohtree->Draw(vardraw.c_str(),mcreal1cut.c_str());

    TH1F *hh_vbf, *hh_vbf_real2, *hh_vbf_real1, *hh_vbf_fake;
    if(hist) {
      hh_vbf = new TH1F("hh_vbf","",nbins,xmin,xmax);
      real2name="hh_vbf_real_"+j.str();
      hh_vbf_real2 = new TH1F(real2name.c_str(),"",nbins,xmin,xmax);
      real1name="hh_vbf_half_"+j.str();
      hh_vbf_real1 = new TH1F(real1name.c_str(),"",nbins,xmin,xmax);
      fakename="hh_vbf_fake_"+j.str();
      hh_vbf_fake = new TH1F(fakename.c_str(),"",nbins,xmin,xmax);
    }
    else {
      hh_vbf = new TH1F("hh_vbf","",13,massLEdges);
      real2name="hh_vbf_real_"+j.str();
      hh_vbf_real2 = new TH1F(real2name.c_str(),"",13,massLEdges);
      real1name="hh_vbf_half_"+j.str();
      hh_vbf_real1 = new TH1F(real1name.c_str(),"",13,massLEdges);
      fakename="hh_vbf_fake_"+j.str();
      hh_vbf_fake = new TH1F(fakename.c_str(),"",13,massLEdges);
    }
    vardraw = var+">>"+fakename;
    hhtree_vbf->Draw(vardraw.c_str(),mcfakecut.c_str());
    vardraw = var+">>"+real2name;
    hhtree_vbf->Draw(vardraw.c_str(),mcreal2cut.c_str());
    vardraw = var+">>"+real1name;
    hhtree_vbf->Draw(vardraw.c_str(),mcreal1cut.c_str());

    TH1F *ttbar, *ttbar_real2, *ttbar_real1, *ttbar_fake;
    if(hist) {
      ttbar = new TH1F("TTbar","",nbins,xmin,xmax);
      real2name="TTbar_real_"+j.str();
      ttbar_real2 = new TH1F(real2name.c_str(),"",nbins,xmin,xmax);
      real1name="TTbar_half_"+j.str();
      ttbar_real1 = new TH1F(real1name.c_str(),"",nbins,xmin,xmax);
      fakename="TTbar_fake_"+j.str();
      ttbar_fake = new TH1F(fakename.c_str(),"",nbins,xmin,xmax);
    }
    else {
      ttbar = new TH1F("TTbar","",13, massLEdges);
      real2name="TTbar_real_"+j.str();
      ttbar_real2 = new TH1F(real2name.c_str(),"",13,massLEdges);
      real1name="TTbar_half_"+j.str();
      ttbar_real1 = new TH1F(real1name.c_str(),"",13,massLEdges);
      fakename="TTbar_fake_"+j.str();
      ttbar_fake = new TH1F(fakename.c_str(),"",13,massLEdges);
    }
    vardraw = var+">>"+fakename;
    tttree->Draw(vardraw.c_str(),mcfakecut.c_str());
    vardraw = var+">>"+real2name;
    tttree->Draw(vardraw.c_str(),mcreal2cut.c_str());
    vardraw = var+">>"+real1name;
    tttree->Draw(vardraw.c_str(),mcreal1cut.c_str());

    smhh->Add(smhh_real2);
    smhh->Add(smhh_real1);
    smhh->Add(smhh_fake);

    ttbar->Add(ttbar_real2);
    ttbar->Add(ttbar_real1);
    ttbar->Add(ttbar_fake);

    Ztt->Add(Ztt_real2);
    Ztt->Add(Ztt_real1);
    Ztt->Add(Ztt_fake);

    Wjets->Add(Wjets_real2);
    Wjets->Add(Wjets_real1);
    Wjets->Add(Wjets_fake);

    ewk->Add(ewk_real2);
    ewk->Add(ewk_real1);
    ewk->Add(ewk_fake);

    vbfh->Add(vbfh_real2);
    vbfh->Add(vbfh_real1);
    vbfh->Add(vbfh_fake);

    ggh->Add(ggh_real2);
    ggh->Add(ggh_real1);
    ggh->Add(ggh_fake);

    assoh->Add(assoh_real2);
    assoh->Add(assoh_real1);
    assoh->Add(assoh_fake);

    hh_vbf->Add(hh_vbf_real2);
    hh_vbf->Add(hh_vbf_real1);
    hh_vbf->Add(hh_vbf_fake);

    //cout << "2 real taus, sig: " << smhh_real2->Integral() << endl;
    //cout << "1 real taus, sig: " << smhh_real1->Integral() << endl;
    //cout << "fake taus, sig:   " << smhh_fake->Integral() << endl;
    //cout << "real taus, ttb: " << ttbar_real->Integral() << endl;
    //cout << "fake taus, ttb: " << ttbar_fake->Integral() << endl;

    std::string dirname="tautau_"+j.str();
    
    TDirectory* lTD = outDC->mkdir(dirname.c_str());
    outDC->cd(lTD->GetPath());
    ttbar->SetName("data_obs");
    ttbar->SetTitle("data_obs");
    ttbar->Write();
    ttbar->SetName("TT");
    ttbar->SetTitle("TT");
    ttbar->Write();
    smhh->SetName("ggHH");
    smhh->SetTitle("ggHH");
    smhh->Write();    
    Ztt->SetName("ZTT");
    Ztt->SetTitle("ZTT");
    Ztt->Write();
    Wjets->SetName("W");
    Wjets->SetTitle("W");
    Wjets->Write();
    ewk->SetName("VV");
    ewk->SetTitle("VV");
    ewk->Write();
    vbfh->SetName("qqH");
    vbfh->SetTitle("qqH");
    vbfh->Write();
    ggh->SetName("ggH");
    ggh->SetTitle("ggH");
    ggh->Write();
    assoh->SetName("assoH");
    assoh->SetTitle("assoH");
    assoh->Write();
    hh_vbf->SetName("qqHH");
    hh_vbf->SetTitle("qqHH");
    hh_vbf->Write();

  }
  outDC->Close();

}
