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

#include "TList.h"
#include "TTree.h"
#include "TCanvas.h"

#include "../Utils/HttStyles.h"

#endif

// RooFit headers
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooFormulaVar.h"
#include "RooSimultaneous.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooExtendPdf.h"
#include "RooKeysPdf.h"
#include "RooPlot.h"

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


void bbtt_test(std::string var="mt2pileup",int nbins=24, double xmin=0, double xmax=600,std::string xtitle="", std::string ytitle="", double sigscale=10)
{

  //TFile *outDC = new TFile("hh_tt_inputs_2.root","RECREATE");

  SetStyle(); gStyle->SetLineStyleString(11,"20 10");
  TH1::SetDefaultSumw2(1);
 
  std::string dir = "/afs/cern.ch/work/j/jlawhorn/public/ntuples/";
  
  std::stringstream scale; scale << sigscale;
  
  //Cut definitions
  double luminosity = 3000;
  std::stringstream lumi; lumi << luminosity;
  std::string objcut = "(tauCat1==1 && tauCat2==1 && ptTau1>45 && ptTau2>45 && ptB1>30 && ptB2>30 && (bTag1==2||bTag1==3||bTag1==6||bTag1==7) && (bTag2==1||bTag2==3||bTag2==6||bTag2==7))";
  std::string jetcut = objcut+"*(ptTrk1>0 && ptTrk2>0)*(mTT>50&&mTT<130)*(mBB1>80&&mBB1<140)*(ptBB1>150)";//*(mHH>300)";
  //signal region
  std::string mccut = jetcut+"*eventWeight*"+lumi.str();
  std::string vbfcut = jetcut+"*eventWeight*49470*0.0632*"+lumi.str();
  std::string sigcut = jetcut+"*eventWeight*"+lumi.str();
  std::string zjetcut = jetcut+"*eventWeight*(eventType!=3&&eventType!=1)*"+lumi.str();
  std::string wjetcut = jetcut+"*eventWeight*(eventType==3&&eventType!=1)*"+lumi.str();
  std::string ewkcut = jetcut+"*eventWeight*(eventType!=1)*"+lumi.str();
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
  std::string vardraw;
  TTree *smhh = (TTree*)hhtree->CopyTree(sigcut.c_str());

  TTree *ztt = (TTree*)vjettree->CopyTree(zjetcut.c_str());
  TTree *ttbar = (TTree*)tttree->CopyTree(mccut.c_str());
  TTree *wjets = (TTree*)vjettree->CopyTree(wjetcut.c_str());
  TTree *ewk = (TTree*) ewktree->CopyTree(ewkcut.c_str());
  TTree *vbfh = (TTree*) vbfhtree->CopyTree(mccut.c_str());
  TTree *ggh = (TTree*) gfhtree->CopyTree(mccut.c_str());
  TTree *assoh = (TTree*) assohtree->CopyTree(mccut.c_str());
  TTree *qcd = (TTree*) qcdtree->CopyTree(mccut.c_str());

  //---------------------------------------------------------------------------

  RooRealVar drawVar(var.c_str(),var.c_str(),xmin,xmax);
  drawVar.setBins(nbins);
  
  RooDataSet *sigDat = new RooDataSet("ggHH", "ggHH", smhh, RooArgSet(drawVar));
  RooDataSet *zttDat   = new RooDataSet("ZTT", "ZTT", ztt, RooArgSet(drawVar));
  RooDataSet *ttbarDat   = new RooDataSet("TT", "TT", ttbar, RooArgSet(drawVar));
  RooDataSet *wjetsDat   = new RooDataSet("W", "W", wjets, RooArgSet(drawVar));
  RooDataSet *ewkDat   = new RooDataSet("VV", "VV", ewk, RooArgSet(drawVar));
  RooDataSet *vbfhDat   = new RooDataSet("qqH", "qqH", vbfh, RooArgSet(drawVar));
  RooDataSet *gghDat   = new RooDataSet("ggH", "ggH", ggh, RooArgSet(drawVar));
  RooDataSet *assohDat   = new RooDataSet("VH", "VH", assoh, RooArgSet(drawVar));
  RooDataSet *qcdDat   = new RooDataSet("QCD", "QCD", qcd, RooArgSet(drawVar));

  TList *list = new TList;
  list->Add(smhh); list->Add(ztt); list->Add(ttbar);
  list->Add(wjets); list->Add(ewk); list->Add(vbfh); 
  list->Add(ggh); list->Add(assoh); list->Add(qcd);
  TTree *all_of_it = TTree::MergeTrees(list);

  RooDataSet *data_obsDat   = new RooDataSet("data_obs", "data_obs", all_of_it, RooArgSet(drawVar));

  RooKeysPdf sig_est("ggHH_shape", "ggHH_shape", drawVar, *sigDat, RooKeysPdf::MirrorBoth,2);
  RooKeysPdf ztt_est("ZTT_shape", "ZTT_shape", drawVar, *zttDat, RooKeysPdf::MirrorBoth,2);
  RooKeysPdf ttbar_est("TT_shape", "TT_shape", drawVar, *ttbarDat, RooKeysPdf::MirrorBoth,2);
  //RooKeysPdf data_obs_est("data_obs_shape", "data_obs_shape", drawVar, *data_obsDat, RooKeysPdf::MirrorBoth,2);
  RooKeysPdf wjets_est("W_shape", "W_shape", drawVar, *wjetsDat, RooKeysPdf::MirrorBoth,2);
  RooKeysPdf ewk_est("VV_shape", "VV_shape", drawVar, *ewkDat, RooKeysPdf::MirrorBoth,2);
  RooKeysPdf vbfh_est("qqH_shape", "qqH_shape", drawVar, *vbfhDat, RooKeysPdf::MirrorBoth,2);
  RooKeysPdf ggh_est("ggH_shape", "ggH_shape", drawVar, *gghDat, RooKeysPdf::MirrorBoth,2);
  RooKeysPdf assoh_est("VH_shape", "VH_shape", drawVar, *assohDat, RooKeysPdf::MirrorBoth,2);
  RooKeysPdf qcd_est("QCD_shape", "QCD_shape", drawVar, *qcdDat, RooKeysPdf::MirrorBoth,2);

  RooWorkspace *w = new RooWorkspace("tautau", "workspace");
  w->import(sig_est);
  w->import(ztt_est);
  w->import(ttbar_est);
  w->import(*data_obsDat);
  w->import(wjets_est);
  w->import(ewk_est);
  w->import(vbfh_est);
  w->import(ggh_est);
  w->import(assoh_est);
  w->import(qcd_est);
  w->writeToFile("hh_tt_inputs_2.root");

}
