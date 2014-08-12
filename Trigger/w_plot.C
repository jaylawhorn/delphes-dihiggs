#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TVector2.h>
#include <TMath.h>
#include <TGraph.h>
#include <TLine.h>
#include <TChain.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <TPaveText.h>

#include "triggerMenu.h"
#include "MitStyleRemix.hh"

#endif

void w_plot() {

  TCanvas *c1 = MakeCanvas("c1", "c1", 600, 600);
  const char* dataset = "CMS Preliminary, W#rightarrow l#nu, at 14 TeV";

  TString we="(e1pt>0 || e2pt>0)";
  TString wm="(m1pt>0 || m2pt>0)";
  TString wt="(t1pt>0 || t2pt>0)";

  TFile *file = new TFile("/afs/cern.ch/work/j/jlawhorn/w.root");
  TTree *tree = (TTree*)file->Get("Events");
  Trigger::Event data;
  tree->SetBranchAddress("Events",&data);

  Float_t tot=tree->GetEntries();

  vector<Float_t> bw; 
  bw.push_back(125); bw.push_back(180);
  bw.push_back(250); bw.push_back(350);

  cout << "Wenu ----" << endl;

  Float_t channel=float(tree->GetEntries(we));
  Float_t analysis=100*tree->GetEntries("((e1pt>25 && abs(e1eta)<2.4)||(e2pt>25 && abs(e2eta)<2.4))*"+we)/channel;

  cout << "Analysis cut: " << analysis << "%" << endl;

  vector<Float_t> v_single_ele;
  vector<Float_t> v_ele_met;
  vector<Float_t> v_any;

  vector<Float_t> v_single_ele_t;
  vector<Float_t> v_ele_met_t;

  v_single_ele_t.push_back(31.2); v_single_ele_t.push_back(27);
  v_single_ele_t.push_back(25.1); v_single_ele_t.push_back(22.8);

  v_ele_met_t.push_back(102); v_ele_met_t.push_back(95);
  v_ele_met_t.push_back(92); v_ele_met_t.push_back(87);

  v_single_ele.push_back(100*tree->GetEntries("(triggerBits126&16)*"+we)/channel);
  v_single_ele.push_back(100*tree->GetEntries("(triggerBits180&16)*"+we)/channel);
  v_single_ele.push_back(100*tree->GetEntries("(triggerBits250&16)*"+we)/channel);
  v_single_ele.push_back(100*tree->GetEntries("(triggerBits350&16)*"+we)/channel);

  v_ele_met.push_back(100*tree->GetEntries("(triggerBits126&131072)*"+we)/channel);
  v_ele_met.push_back(100*tree->GetEntries("(triggerBits180&131072)*"+we)/channel);
  v_ele_met.push_back(100*tree->GetEntries("(triggerBits250&131072)*"+we)/channel);
  v_ele_met.push_back(100*tree->GetEntries("(triggerBits350&131072)*"+we)/channel);
  
  v_any.push_back(100*tree->GetEntries("((triggerBits126&16)||(triggerBits126&131072))*"+we)/channel);
  v_any.push_back(100*tree->GetEntries("((triggerBits180&16)||(triggerBits180&131072))*"+we)/channel);
  v_any.push_back(100*tree->GetEntries("((triggerBits250&16)||(triggerBits250&131072))*"+we)/channel);
  v_any.push_back(100*tree->GetEntries("((triggerBits350&16)||(triggerBits350&131072))*"+we)/channel);

  TGraph *single_ele = new TGraph(4, &bw[0], &v_single_ele[0]);
  InitGraph(single_ele, "Total Bandwith [kHz]", "Acceptance [%]",kRed);
  single_ele->Draw("ap");
  CMSPrelim(dataset, "W#rightarrow e#nu", 0.16, 0.835);
  c1->SaveAs("wenu_single_ele.png");

  TGraph *single_ele2 = new TGraph(4, &v_single_ele_t[0], &v_single_ele[0]);
  InitGraph(single_ele2, "Threshold [GeV]", "Acceptance [%]",kRed);
  single_ele2->Draw("ap");
  CMSPrelim(dataset, "W#rightarrow e#nu", 0.16, 0.835);
  c1->SaveAs("wenu_single_ele2.png");

  TGraph *ele_met = new TGraph(4, &bw[0], &v_ele_met[0]);
  InitGraph(ele_met, "Total Bandwidth [kHz]", "Acceptance [%]", kRed);
  ele_met->Draw("ap");
  CMSPrelim(dataset, "W#rightarrow e#nu", 0.16, 0.835);
  c1->SaveAs("wenu_ele_met.png");

  TGraph *ele_met2 = new TGraph(4, &v_ele_met_t[0], &v_ele_met[0]);
  InitGraph(ele_met2, "MET Threshold [GeV]", "Acceptance [%]", kRed);
  ele_met2->Draw("ap");
  CMSPrelim(dataset, "W#rightarrow e#nu", 0.16, 0.835);
  c1->SaveAs("wenu_ele_met2.png");

  TGraph *either =new TGraph(4, &bw[0], &v_any[0]);
  InitGraph(either, "Total Bandwith [kHz]", "Acceptance [%]",kBlue);
  either->Draw("ap");
  CMSPrelim(dataset, "W#rightarrow e#nu", 0.16, 0.835);
  c1->SaveAs("wenu_inc.png");

  v_single_ele.clear();
  v_ele_met.clear();
  v_any.clear();
  delete(single_ele);
  delete(ele_met);
  delete(either);

  cout << "Wmunu ---" << endl;

  channel=float(tree->GetEntries(wm));
  analysis=100*tree->GetEntries("((m1pt>25 && abs(m1eta)<2.4)||(m2pt>25 && abs(m2eta)<2.4))*"+wm)/channel;
  
  cout << "Analysis cut: " << analysis << "%" << endl;

  vector<Float_t> v_single_mu;
  vector<Float_t> v_mu_met;
  vector<Float_t> v_single_mu_t;
  vector<Float_t> v_mu_met_t;

  v_single_mu_t.push_back(21.1); v_single_mu_t.push_back(18);
  v_single_mu_t.push_back(17.3); v_single_mu_t.push_back(15.5);

  v_mu_met_t.push_back(68.6); v_mu_met_t.push_back(95);
  v_mu_met_t.push_back(61); v_mu_met_t.push_back(57);

  v_single_mu.push_back(100*tree->GetEntries("(triggerBits126&1)*"+wm)/channel);
  v_single_mu.push_back(100*tree->GetEntries("(triggerBits180&1)*"+wm)/channel);
  v_single_mu.push_back(100*tree->GetEntries("(triggerBits250&1)*"+wm)/channel);
  v_single_mu.push_back(100*tree->GetEntries("(triggerBits350&1)*"+wm)/channel);

  v_mu_met.push_back(100*tree->GetEntries("(triggerBits126&262144)*"+wm)/channel);
  v_mu_met.push_back(100*tree->GetEntries("(triggerBits180&262144)*"+wm)/channel);
  v_mu_met.push_back(100*tree->GetEntries("(triggerBits250&262144)*"+wm)/channel);
  v_mu_met.push_back(100*tree->GetEntries("(triggerBits350&262144)*"+wm)/channel);
  
  v_any.push_back(100*tree->GetEntries("((triggerBits126&1)||(triggerBits126&262144))*"+wm)/channel);
  v_any.push_back(100*tree->GetEntries("((triggerBits180&1)||(triggerBits180&262144))*"+wm)/channel);
  v_any.push_back(100*tree->GetEntries("((triggerBits250&1)||(triggerBits250&262144))*"+wm)/channel);
  v_any.push_back(100*tree->GetEntries("((triggerBits350&1)||(triggerBits350&262144))*"+wm)/channel);

  TGraph *single_mu = new TGraph(4, &bw[0], &v_single_mu[0]);
  InitGraph(single_mu, "Total Bandwith [kHz]", "Acceptance [%]",kRed);
  single_mu->Draw("ap");
  CMSPrelim(dataset, "W#rightarrow #mu#nu", 0.16, 0.835);
  c1->SaveAs("wmnu_single_mu.png");

  TGraph *single_mu2 = new TGraph(4, &v_single_mu_t[0], &v_single_mu[0]);
  InitGraph(single_mu2, "Threshold [GeV]", "Acceptance [%]",kRed);
  single_mu2->Draw("ap");
  CMSPrelim(dataset, "W#rightarrow #mu#nu", 0.16, 0.835);
  c1->SaveAs("wmnu_single_mu2.png");

  TGraph *mu_met = new TGraph(4, &bw[0], &v_mu_met[0]);
  InitGraph(mu_met, "Total Bandwidth [kHz]", "Acceptance [%]", kRed);
  mu_met->Draw("ap");
  CMSPrelim(dataset, "W#rightarrow #mu#nu", 0.16, 0.835);
  c1->SaveAs("wmnu_mu_met.png");

  TGraph *mu_met2 = new TGraph(4, &v_mu_met_t[0], &v_mu_met[0]);
  InitGraph(mu_met2, "MET Threshold [GeV]", "Acceptance [%]", kRed);
  mu_met2->Draw("ap");
  CMSPrelim(dataset, "W#rightarrow #mu#nu", 0.16, 0.835);
  c1->SaveAs("wmnu_mu_met2.png");

  either = new TGraph(4, &bw[0], &v_any[0]);
  InitGraph(either, "Total Bandwith [kHz]", "Acceptance [%]",kBlue);
  either->Draw("ap");
  CMSPrelim(dataset, "W#rightarrow #mu#nu", 0.16, 0.835);
  c1->SaveAs("wmnu_inc.png");

}
