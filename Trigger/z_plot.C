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

void z_plot() {

  TCanvas *c1 = MakeCanvas("c1", "c1", 600, 600);
  const char* dataset = "CMS Preliminary, Z#rightarrow ll, at 14 TeV";

  TString zee="(e1pt>0 && e2pt>0)";
  TString zmm="(m1pt>0 && m2pt>0)";
  TString ztt="(t1pt>0 && t2pt>0)";

  TFile *file = new TFile("/afs/cern.ch/work/j/jlawhorn/z.root");
  TTree *tree = (TTree*)file->Get("Events");
  Trigger::Event data;
  tree->SetBranchAddress("Events",&data);

  Float_t tot=tree->GetEntries();

  vector<Float_t> bw; 
  bw.push_back(125); bw.push_back(180);
  bw.push_back(250); bw.push_back(350);

  cout << "Zee ----" << endl;

  Float_t channel=float(tree->GetEntries(zee));
  Float_t analysis=100*tree->GetEntries("((e1pt>25 && abs(e1eta)<2.4)&&(e2pt>25 && abs(e2eta)<2.4))*"+zee)/channel;

  cout << "Analysis cut: " << analysis << "%" << endl;

  vector<Float_t> v_single_ele;
  vector<Float_t> v_double_ele;
  vector<Float_t> v_single_ele_t;
  vector<Float_t> v_double_ele_t;
  vector<Float_t> v_any;

  v_single_ele_t.push_back(31.2); v_single_ele_t.push_back(27);
  v_single_ele_t.push_back(25.1); v_single_ele_t.push_back(22.8);

  v_double_ele_t.push_back(23.3); v_double_ele_t.push_back(22);
  v_double_ele_t.push_back(21.5); v_double_ele_t.push_back(20.6);

  v_single_ele.push_back(100*tree->GetEntries("(triggerBits126&16)*"+zee)/channel);
  v_single_ele.push_back(100*tree->GetEntries("(triggerBits180&16)*"+zee)/channel);
  v_single_ele.push_back(100*tree->GetEntries("(triggerBits250&16)*"+zee)/channel);
  v_single_ele.push_back(100*tree->GetEntries("(triggerBits350&16)*"+zee)/channel);

  v_double_ele.push_back(100*tree->GetEntries("(triggerBits126&64)*"+zee)/channel);
  v_double_ele.push_back(100*tree->GetEntries("(triggerBits180&64)*"+zee)/channel);
  v_double_ele.push_back(100*tree->GetEntries("(triggerBits250&64)*"+zee)/channel);
  v_double_ele.push_back(100*tree->GetEntries("(triggerBits350&64)*"+zee)/channel);

  v_any.push_back(100*tree->GetEntries("((triggerBits126&16)||(triggerBits126&64))*"+zee)/channel);
  v_any.push_back(100*tree->GetEntries("((triggerBits180&16)||(triggerBits180&64))*"+zee)/channel);
  v_any.push_back(100*tree->GetEntries("((triggerBits250&16)||(triggerBits250&64))*"+zee)/channel);
  v_any.push_back(100*tree->GetEntries("((triggerBits350&16)||(triggerBits350&64))*"+zee)/channel);

  TGraph *single_ele = new TGraph(4, &bw[0], &v_single_ele[0]);
  InitGraph(single_ele, "Total Bandwith [kHz]", "Acceptance [%]",kRed);
  single_ele->Draw("ap");
  CMSPrelim(dataset, "Z#rightarrow ee", 0.16, 0.835);
  c1->SaveAs("zee_single_ele.png");

  TGraph *single_ele2 = new TGraph(4, &v_single_ele_t[0], &v_single_ele[0]);
  InitGraph(single_ele2, "Threshold [GeV]", "Acceptance [%]",kRed);
  single_ele2->Draw("ap");
  CMSPrelim(dataset, "Z#rightarrow ee", 0.16, 0.835);
  c1->SaveAs("zee_single_ele2.png");

  TGraph *double_ele = new TGraph(4, &bw[0], &v_double_ele[0]);
  InitGraph(double_ele, "Total Bandwidth [kHz]", "Acceptance [%]", kRed);
  double_ele->Draw("ap");
  CMSPrelim(dataset, "Z#rightarrow ee", 0.16, 0.835);
  c1->SaveAs("zee_double_ele.png");

  TGraph *double_ele2 = new TGraph(4, &v_double_ele_t[0], &v_double_ele[0]);
  InitGraph(double_ele2, "Leading e Threshold [GeV]", "Acceptance [%]", kRed);
  double_ele2->Draw("ap");
  CMSPrelim(dataset, "Z#rightarrow ee", 0.16, 0.835);
  c1->SaveAs("zee_double_ele2.png");

  TGraph *either =new TGraph(4, &bw[0], &v_any[0]);
  InitGraph(either, "Total Bandwith [kHz]", "Acceptance [%]",kBlue);
  either->Draw("ap");
  CMSPrelim(dataset, "Z#rightarrow ee", 0.16, 0.835);
  c1->SaveAs("zee_inc.png");

  v_any.clear();
  delete(either);

  cout << "Zmm ---" << endl;

  channel=float(tree->GetEntries(zmm));
  analysis=100*tree->GetEntries("((m1pt>25 && abs(m1eta)<2.4)&&(m2pt>25 && abs(m2eta)<2.4))*"+zmm)/channel;
  
  cout << "Analysis cut: " << analysis << "%" << endl;

  vector<Float_t> v_single_mu;
  vector<Float_t> v_double_mu;
  vector<Float_t> v_single_mu_t;
  vector<Float_t> v_double_mu_t;

  v_single_mu_t.push_back(21.1); v_single_mu_t.push_back(18);
  v_single_mu_t.push_back(17.3); v_single_mu_t.push_back(15.5);

  v_double_mu_t.push_back(18.4); v_double_mu_t.push_back(14);
  v_double_mu_t.push_back(11.9); v_double_mu_t.push_back(10.4);

  v_single_mu.push_back(100*tree->GetEntries("(triggerBits126&1)*"+zmm)/channel);
  v_single_mu.push_back(100*tree->GetEntries("(triggerBits180&1)*"+zmm)/channel);
  v_single_mu.push_back(100*tree->GetEntries("(triggerBits250&1)*"+zmm)/channel);
  v_single_mu.push_back(100*tree->GetEntries("(triggerBits350&1)*"+zmm)/channel);

  v_double_mu.push_back(100*tree->GetEntries("(triggerBits126&2)*"+zmm)/channel);
  v_double_mu.push_back(100*tree->GetEntries("(triggerBits180&2)*"+zmm)/channel);
  v_double_mu.push_back(100*tree->GetEntries("(triggerBits250&2)*"+zmm)/channel);
  v_double_mu.push_back(100*tree->GetEntries("(triggerBits350&2)*"+zmm)/channel);

  v_any.push_back(100*tree->GetEntries("((triggerBits126&1)||(triggerBits126&2))*"+zmm)/channel);
  v_any.push_back(100*tree->GetEntries("((triggerBits180&1)||(triggerBits180&2))*"+zmm)/channel);
  v_any.push_back(100*tree->GetEntries("((triggerBits250&1)||(triggerBits250&2))*"+zmm)/channel);
  v_any.push_back(100*tree->GetEntries("((triggerBits350&1)||(triggerBits350&2))*"+zmm)/channel);

  TGraph *single_mu = new TGraph(4, &bw[0], &v_single_mu[0]);
  InitGraph(single_mu, "Total Bandwith [kHz]", "Acceptance [%]",kRed);
  single_mu->Draw("ap");
  CMSPrelim(dataset, "Z#rightarrow #mu#mu", 0.16, 0.835);
  c1->SaveAs("zmm_single_mu.png");

  TGraph *single_mu2 = new TGraph(4, &v_single_mu_t[0], &v_single_mu[0]);
  InitGraph(single_mu2, "Threshold [GeV]", "Acceptance [%]",kRed);
  single_mu2->Draw("ap");
  CMSPrelim(dataset, "Z#rightarrow #mu#mu", 0.16, 0.835);
  c1->SaveAs("zmm_single_mu2.png");

  TGraph *double_mu = new TGraph(4, &bw[0], &v_double_mu[0]);
  InitGraph(double_mu, "Total Bandwith [kHz]", "Acceptance [%]",kRed);
  double_mu->Draw("ap");
  CMSPrelim(dataset, "Z#rightarrow #mu#mu", 0.16, 0.835);
  c1->SaveAs("zmm_double_mu.png");

  TGraph *double_mu2 = new TGraph(4, &v_double_mu_t[0], &v_double_mu[0]);
  InitGraph(double_mu2, "Leading #mu Threshold [GeV]", "Acceptance [%]",kRed);
  double_mu2->Draw("ap");
  CMSPrelim(dataset, "Z#rightarrow #mu#mu", 0.16, 0.835);
  c1->SaveAs("zmm_double_mu2.png");

  either = new TGraph(4, &bw[0], &v_any[0]);
  InitGraph(either, "Total Bandwith [kHz]", "Acceptance [%]",kBlue);
  either->Draw("ap");
  CMSPrelim(dataset, "Z#rightarrow #mu#mu", 0.16, 0.835);
  c1->SaveAs("zmm_inc.png");

}
