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

void susy_plot() {

  TCanvas *c1 = MakeCanvas("c1", "c1", 600, 600);
  const char* dataset = "CMS Preliminary, SUSY NM3, at 14 TeV";

  TFile *file = new TFile("/afs/cern.ch/work/j/jlawhorn/susy.root");
  TTree *tree = (TTree*)file->Get("Events");
  Trigger::Event data;
  tree->SetBranchAddress("Events",&data);

  Float_t tot=tree->GetEntries();

  vector<Float_t> bw; 
  bw.push_back(125); bw.push_back(180);
  bw.push_back(250); bw.push_back(350);

  cout << 100*tree->GetEntries("(triggerBits350&1)")/tot << endl;
  cout << 100*tree->GetEntries("(triggerBits350&2)")/tot << endl;
  cout << 100*tree->GetEntries("(triggerBits350&4)")/tot << endl;
  cout << 100*tree->GetEntries("(triggerBits350&8)")/tot << endl;
  cout << 100*tree->GetEntries("(triggerBits350&16)")/tot << endl;
  cout << 100*tree->GetEntries("(triggerBits350&32)")/tot << endl;
  cout << 100*tree->GetEntries("(triggerBits350&64)")/tot << endl;
  cout << 100*tree->GetEntries("(triggerBits350&128)")/tot << endl;
  cout << 100*tree->GetEntries("(triggerBits350&256)")/tot << endl;
  cout << 100*tree->GetEntries("(triggerBits350&512)")/tot << endl;
  cout << 100*tree->GetEntries("(triggerBits350&1024)")/tot << endl;
  cout << 100*tree->GetEntries("(triggerBits350&2048)")/tot << endl;
  cout << 100*tree->GetEntries("(triggerBits350&4096)")/tot << endl;
  cout << 100*tree->GetEntries("(triggerBits350&8192)")/tot << endl;
  cout << 100*tree->GetEntries("(triggerBits350&16384)")/tot << endl;
  cout << 100*tree->GetEntries("(triggerBits350&32768)")/tot << endl;
  cout << 100*tree->GetEntries("(triggerBits350&65536)")/tot << endl;
  cout << 100*tree->GetEntries("(triggerBits350&131072)")/tot << endl;
  cout << 100*tree->GetEntries("(triggerBits350&262114)")/tot << endl;
  cout << 100*tree->GetEntries("(triggerBits350&524288)")/tot << endl;

  vector<Float_t> v_single_ele;
  vector<Float_t> v_single_jet;
  vector<Float_t> v_single_mu;
  vector<Float_t> v_mu_met;
  vector<Float_t> v_ele_met;
  vector<Float_t> v_ht;
  vector<Float_t> v_single_ele_t;
  vector<Float_t> v_single_jet_t;
  vector<Float_t> v_single_mu_t;
  vector<Float_t> v_mu_met_t;
  vector<Float_t> v_ele_met_t;
  vector<Float_t> v_ht_t;
  vector<Float_t> v_any;

  v_single_ele_t.push_back(31.2); v_single_ele_t.push_back(27);
  v_single_ele_t.push_back(25.1); v_single_ele_t.push_back(22.8);

  v_single_jet_t.push_back(188); v_single_jet_t.push_back(173);
  v_single_jet_t.push_back(167); v_single_jet_t.push_back(158);

  v_single_mu_t.push_back(21.1); v_single_mu_t.push_back(18);
  v_single_mu_t.push_back(17.3); v_single_mu_t.push_back(15.5);

  v_mu_met_t.push_back(68.6); v_mu_met_t.push_back(95);
  v_mu_met_t.push_back(61); v_mu_met_t.push_back(57);

  v_ele_met_t.push_back(102); v_ele_met_t.push_back(95);
  v_ele_met_t.push_back(92); v_ele_met_t.push_back(87);

  v_ht_t.push_back(384); v_ht_t.push_back(350);
  v_ht_t.push_back(340); v_ht_t.push_back(317);

  v_single_ele.push_back(100*tree->GetEntries("(triggerBits126&16)")/tot);
  v_single_ele.push_back(100*tree->GetEntries("(triggerBits180&16)")/tot);
  v_single_ele.push_back(100*tree->GetEntries("(triggerBits250&16)")/tot);
  v_single_ele.push_back(100*tree->GetEntries("(triggerBits350&16)")/tot);

  v_single_mu.push_back(100*tree->GetEntries("(triggerBits126&1)")/tot);
  v_single_mu.push_back(100*tree->GetEntries("(triggerBits180&1)")/tot);
  v_single_mu.push_back(100*tree->GetEntries("(triggerBits250&1)")/tot);
  v_single_mu.push_back(100*tree->GetEntries("(triggerBits350&1)")/tot);

  v_single_jet.push_back(100*tree->GetEntries("(triggerBits126&4096)")/tot);
  v_single_jet.push_back(100*tree->GetEntries("(triggerBits180&4096)")/tot);
  v_single_jet.push_back(100*tree->GetEntries("(triggerBits250&4096)")/tot);
  v_single_jet.push_back(100*tree->GetEntries("(triggerBits350&4096)")/tot);

  v_ele_met.push_back(100*tree->GetEntries("(triggerBits126&131072)")/tot);
  v_ele_met.push_back(100*tree->GetEntries("(triggerBits180&131072)")/tot);
  v_ele_met.push_back(100*tree->GetEntries("(triggerBits250&131072)")/tot);
  v_ele_met.push_back(100*tree->GetEntries("(triggerBits350&131072)")/tot);

  v_mu_met.push_back(100*tree->GetEntries("(triggerBits126&262144)")/tot);
  v_mu_met.push_back(100*tree->GetEntries("(triggerBits180&262144)")/tot);
  v_mu_met.push_back(100*tree->GetEntries("(triggerBits250&262144)")/tot);
  v_mu_met.push_back(100*tree->GetEntries("(triggerBits350&262144)")/tot);

  v_ht.push_back(100*tree->GetEntries("(triggerBits126&524288)")/tot);
  v_ht.push_back(100*tree->GetEntries("(triggerBits180&524288)")/tot);
  v_ht.push_back(100*tree->GetEntries("(triggerBits250&524288)")/tot);
  v_ht.push_back(100*tree->GetEntries("(triggerBits350&524288)")/tot);
  
  v_any.push_back(100*tree->GetEntries("((triggerBits126&16)||(triggerBits126&1)||(triggerBits126&4096)||(triggerBits126&131072)||(triggerBits126&262144)||(triggerBits126&524288))")/tot);
  v_any.push_back(100*tree->GetEntries("((triggerBits180&16)||(triggerBits180&1)||(triggerBits180&4096)||(triggerBits180&131072)||(triggerBits180&262144)||(triggerBits180&524288))")/tot);
  v_any.push_back(100*tree->GetEntries("((triggerBits250&16)||(triggerBits250&1)||(triggerBits250&4096)||(triggerBits250&131072)||(triggerBits250&262144)||(triggerBits250&524288))")/tot);
  v_any.push_back(100*tree->GetEntries("((triggerBits350&16)||(triggerBits350&1)||(triggerBits350&4096)||(triggerBits350&131072)||(triggerBits350&262144)||(triggerBits350&524288))")/tot);

  TGraph *single_ele = new TGraph(4, &bw[0], &v_single_ele[0]);
  InitGraph(single_ele, "Total Bandwith [kHz]", "Acceptance [%]",kRed);
  single_ele->Draw("ap");
  CMSPrelim(dataset, "", 0.16, 0.835);
  c1->SaveAs("susy_single_ele.png");

  TGraph *single_ele2 = new TGraph(4, &v_single_ele_t[0], &v_single_ele[0]);
  InitGraph(single_ele2, "Threshold [GeV]", "Acceptance [%]",kRed);
  single_ele2->Draw("ap");
  CMSPrelim(dataset, "", 0.16, 0.835);
  c1->SaveAs("susy_single_ele2.png");

  TGraph *single_mu = new TGraph(4, &bw[0], &v_single_mu[0]);
  InitGraph(single_mu, "Total Bandwith [kHz]", "Acceptance [%]",kRed);
  single_mu->Draw("ap");
  CMSPrelim(dataset, "", 0.16, 0.835);
  c1->SaveAs("susy_single_mu.png");

  TGraph *single_mu2 = new TGraph(4, &v_single_mu_t[0], &v_single_mu[0]);
  InitGraph(single_mu2, "Threshold [GeV]", "Acceptance [%]",kRed);
  single_mu2->Draw("ap");
  CMSPrelim(dataset, "", 0.16, 0.835);
  c1->SaveAs("susy_single_mu2.png");

  TGraph *single_jet = new TGraph(4, &bw[0], &v_single_jet[0]);
  InitGraph(single_jet, "Total Bandwith [kHz]", "Acceptance [%]",kRed);
  single_jet->Draw("ap");
  CMSPrelim(dataset, "", 0.16, 0.835);
  c1->SaveAs("susy_single_jet.png");

  TGraph *single_jet2 = new TGraph(4, &v_single_jet_t[0], &v_single_jet[0]);
  InitGraph(single_jet2, "Threshold [GeV]", "Acceptance [%]",kRed);
  single_jet2->Draw("ap");
  CMSPrelim(dataset, "", 0.16, 0.835);
  c1->SaveAs("susy_single_jet2.png");

  TGraph *ele_met = new TGraph(4, &bw[0], &v_ele_met[0]);
  InitGraph(ele_met, "Total Bandwidth [kHz]", "Acceptance [%]", kRed);
  ele_met->Draw("ap");
  CMSPrelim(dataset, "", 0.16, 0.835);
  c1->SaveAs("susy_ele_met.png");

  TGraph *ele_met2 = new TGraph(4, &v_ele_met_t[0], &v_ele_met[0]);
  InitGraph(ele_met2, "MET Threshold [GeV]", "Acceptance [%]", kRed);
  ele_met2->Draw("ap");
  CMSPrelim(dataset, "", 0.16, 0.835);
  c1->SaveAs("susy_ele_met2.png");

  TGraph *mu_met = new TGraph(4, &bw[0], &v_mu_met[0]);
  InitGraph(mu_met, "Total Bandwidth [kHz]", "Acceptance [%]", kRed);
  mu_met->Draw("ap");
  CMSPrelim(dataset, "", 0.16, 0.835);
  c1->SaveAs("susy_mu_met.png");

  TGraph *mu_met2 = new TGraph(4, &v_mu_met_t[0], &v_mu_met[0]);
  InitGraph(mu_met2, "MET Threshold [GeV]", "Acceptance [%]", kRed);
  mu_met2->Draw("ap");
  CMSPrelim(dataset, "", 0.16, 0.835);
  c1->SaveAs("susy_mu_met2.png");

  TGraph *ht = new TGraph(4, &bw[0], &v_ht[0]);
  InitGraph(ht, "Total Bandwith [kHz]", "Acceptance [%]",kRed);
  ht->Draw("ap");
  CMSPrelim(dataset, "", 0.16, 0.835);
  c1->SaveAs("susy_ht.png");

  TGraph *ht2 = new TGraph(4, &v_ht_t[0], &v_ht[0]);
  InitGraph(ht2, "Threshold [GeV]", "Acceptance [%]",kRed);
  ht2->Draw("ap");
  CMSPrelim(dataset, "", 0.16, 0.835);
  c1->SaveAs("susy_ht2.png");

  TGraph *either =new TGraph(4, &bw[0], &v_any[0]);
  InitGraph(either, "Total Bandwith [kHz]", "Acceptance [%]",kBlue);
  either->Draw("ap");
  CMSPrelim(dataset, "", 0.16, 0.835);
  c1->SaveAs("susy_inc.png");


}
