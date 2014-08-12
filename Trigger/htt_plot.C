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

void htt_plot() {

  TCanvas *c1 = MakeCanvas("c1", "c1", 600, 600);
  const char* dataset = "CMS Preliminary, H#rightarrow #tau#tau, at 14 TeV";

  TFile *file = new TFile("/afs/cern.ch/work/j/jlawhorn/htt.root");
  TTree *tree = (TTree*)file->Get("Events");
  Trigger::Event data;
  tree->SetBranchAddress("Events",&data);

  Float_t tot=tree->GetEntries();

  vector<Float_t> bw; 
  bw.push_back(125); bw.push_back(180);
  bw.push_back(250); bw.push_back(350);

  TString tt="t1pt>0 && t2pt>0";
  TString mt="m1pt>0 && (t1pt>0 || t2pt>0)";
  TString et="e1pt>0 && (t1pt>0 || t2pt>0)";
  TString em="m1pt>0 && e1pt>0";
  TString ee="e1pt>0 && e2pt>0";
  TString mm="m1pt>0 && m2pt>0";
  /*  
  cout << "Total entries: " << tot << endl;
  cout << "Tau-tau:       " << 100*tree->GetEntries(tt)/tot << "%" << endl;
  cout << "Mu-tau:        " << 100*tree->GetEntries(mt)/tot << "%" << endl;
  cout << "Ele-tau:       " << 100*tree->GetEntries(et)/tot << "%" << endl;
  cout << "Ele-mu:        " << 100*tree->GetEntries(em)/tot << "%" << endl;
  cout << "ee:            " << 100*tree->GetEntries(ee)/tot << "%" << endl;
  cout << "mm:            " << 100*tree->GetEntries(mm)/tot << "%" << endl;
  cout << "Nothing:        " << 100*tree->GetEntries("!("+tt+")*!("+mt+")*!("+et+")*!("+em+")*!("+ee+")*!("+mm+")")/tot << "%" << endl;
  cout << endl;
  */
  cout << "tt ----" << endl;

  Float_t channel=float(tree->GetEntries(tt));
  Float_t analysis=100*tree->GetEntries("(t1pt>45 && t2pt>45 && abs(t1eta)<2.4 && abs(t2eta)<2.4)*"+tt)/channel;

  cout << "Analysis: " << analysis << "%" << endl;
  vector<Float_t> v_single_tau;
  vector<Float_t> v_single_tau_t;
  v_single_tau_t.push_back(101); v_single_tau_t.push_back(88);
  v_single_tau_t.push_back(83.9); v_single_tau_t.push_back(76);
  vector<Float_t> v_double_tau;
  vector<Float_t> v_double_tau_t;
  v_double_tau_t.push_back(63.8); v_double_tau_t.push_back(56);
  v_double_tau_t.push_back(54); v_double_tau_t.push_back(49.9);
  vector<Float_t> v_either;

  v_single_tau.push_back(100*tree->GetEntries("(triggerBits126&256)*"+tt)/channel);
  v_single_tau.push_back(100*tree->GetEntries("(triggerBits180&256)*"+tt)/channel);
  v_single_tau.push_back(100*tree->GetEntries("(triggerBits250&256)*"+tt)/channel);
  v_single_tau.push_back(100*tree->GetEntries("(triggerBits350&256)*"+tt)/channel);

  v_double_tau.push_back(100*tree->GetEntries("(triggerBits126&512)*"+tt)/channel);
  v_double_tau.push_back(100*tree->GetEntries("(triggerBits180&512)*"+tt)/channel);
  v_double_tau.push_back(100*tree->GetEntries("(triggerBits250&512)*"+tt)/channel);
  v_double_tau.push_back(100*tree->GetEntries("(triggerBits350&512)*"+tt)/channel);
  
  v_either.push_back(100*tree->GetEntries("((triggerBits126&512)||(triggerBits126&256))*"+tt)/channel);
  v_either.push_back(100*tree->GetEntries("((triggerBits180&512)||(triggerBits180&256))*"+tt)/channel);
  v_either.push_back(100*tree->GetEntries("((triggerBits250&512)||(triggerBits250&256))*"+tt)/channel);
  v_either.push_back(100*tree->GetEntries("((triggerBits350&512)||(triggerBits350&256))*"+tt)/channel);

  TGraph *single_tau = new TGraph(4, &bw[0], &v_single_tau[0]);
  InitGraph(single_tau, "Total Bandwith [kHz]", "Acceptance [%]",kRed);
  single_tau->Draw("ap");
  CMSPrelim(dataset, "#tau_{h}#tau_{h}", 0.16, 0.835);
  c1->SaveAs("thth_single_tau.png");

  TGraph *single_tau2 = new TGraph(4, &v_single_tau_t[0], &v_single_tau[0]);
  InitGraph(single_tau2, "Threshold [GeV]", "Acceptance [%]",kRed);
  single_tau2->Draw("ap");
  CMSPrelim(dataset, "#tau_{h}#tau_{h}", 0.16, 0.835);
  c1->SaveAs("thth_single_tau2.png");

  TGraph *double_tau = new TGraph(4, &bw[0], &v_double_tau[0]);
  InitGraph(double_tau, "Total Bandwith [kHz]", "Acceptance [%]",kRed);
  double_tau->Draw("ap");
  CMSPrelim(dataset, "#tau_{h}#tau_{h}", 0.16, 0.835);
  c1->SaveAs("thth_double_tau.png");

  TGraph *double_tau2 = new TGraph(4, &v_double_tau_t[0], &v_double_tau[0]);
  InitGraph(double_tau2, "Leading #tau Threshold [GeV]", "Acceptance [%]",kRed);
  double_tau2->Draw("ap");
  CMSPrelim(dataset, "#tau_{h}#tau_{h}", 0.16, 0.835);
  c1->SaveAs("thth_double_tau2.png");

  TGraph *either = new TGraph(4, &bw[0], &v_either[0]);
  InitGraph(either, "Total Bandwith [kHz]", "Acceptance [%]",kBlue);
  either->Draw("ap");
  CMSPrelim(dataset, "#tau_{h}#tau_{h}", 0.16, 0.835);
  c1->SaveAs("thth_inc.png");

  v_single_tau.clear();
  v_double_tau.clear();
  v_either.clear();
  delete(single_tau);
  delete(single_tau2);
  delete(double_tau);
  delete(double_tau2);
  delete(either);

  cout << "mt ----" << endl;

  channel=float(tree->GetEntries(mt));
  analysis=100*tree->GetEntries("(((t1pt>30 && abs(t1eta)<2.4)||(t2pt>30 && abs(t2eta)<2.4))&&((m1pt>20 && abs(m1eta)<2.4)||(m2pt>20 && abs(m2eta)<2.4)))*"+mt)/channel;
  cout << "Analysis: " << analysis << "%" << endl;

  vector<Float_t> v_single_mu;
  vector<Float_t> v_single_mu_t;
  v_single_mu_t.push_back(21.1); v_single_mu_t.push_back(18);
  v_single_mu_t.push_back(17.3); v_single_mu_t.push_back(15.5);
  vector<Float_t> v_mu_tau;
  vector<Float_t> v_mu_tau_t;
  v_mu_tau_t.push_back(51.9); v_mu_tau_t.push_back(45);
  v_mu_tau_t.push_back(43.4); v_mu_tau_t.push_back(40.8);

  v_single_tau.push_back(100*tree->GetEntries("(triggerBits126&256)*"+mt)/channel);
  v_single_tau.push_back(100*tree->GetEntries("(triggerBits180&256)*"+mt)/channel);
  v_single_tau.push_back(100*tree->GetEntries("(triggerBits250&256)*"+mt)/channel);
  v_single_tau.push_back(100*tree->GetEntries("(triggerBits350&256)*"+mt)/channel);

  v_single_mu.push_back(100*tree->GetEntries("(triggerBits126&1)*"+mt)/channel);
  v_single_mu.push_back(100*tree->GetEntries("(triggerBits180&1)*"+mt)/channel);
  v_single_mu.push_back(100*tree->GetEntries("(triggerBits250&1)*"+mt)/channel);
  v_single_mu.push_back(100*tree->GetEntries("(triggerBits350&1)*"+mt)/channel);

  v_mu_tau.push_back(100*tree->GetEntries("(triggerBits126&2048)*"+mt)/channel);
  v_mu_tau.push_back(100*tree->GetEntries("(triggerBits180&2048)*"+mt)/channel);
  v_mu_tau.push_back(100*tree->GetEntries("(triggerBits250&2048)*"+mt)/channel);
  v_mu_tau.push_back(100*tree->GetEntries("(triggerBits350&2048)*"+mt)/channel);

  v_either.push_back(100*tree->GetEntries("((triggerBits126&256)||(triggerBits126&1)||(triggerBits126&2048))*"+mt)/channel);
  v_either.push_back(100*tree->GetEntries("((triggerBits180&256)||(triggerBits180&1)||(triggerBits180&2048))*"+mt)/channel);
  v_either.push_back(100*tree->GetEntries("((triggerBits250&256)||(triggerBits250&1)||(triggerBits250&2048))*"+mt)/channel);
  v_either.push_back(100*tree->GetEntries("((triggerBits350&256)||(triggerBits350&1)||(triggerBits350&2048))*"+mt)/channel);

  single_tau = new TGraph(4, &bw[0], &v_single_tau[0]);
  InitGraph(single_tau, "Total Bandwith [kHz]", "Acceptance [%]",kRed);
  single_tau->Draw("ap");
  CMSPrelim(dataset, "#mu#tau_{h}", 0.16, 0.835);
  c1->SaveAs("thtm_single_tau.png");

  single_tau2 = new TGraph(4, &v_single_tau_t[0], &v_single_tau[0]);
  InitGraph(single_tau2, "Threshold [GeV]", "Acceptance [%]",kRed);
  single_tau2->Draw("ap");
  CMSPrelim(dataset, "#mu#tau_{h}", 0.16, 0.835);
  c1->SaveAs("thtm_single_tau2.png");

  TGraph *single_mu = new TGraph(4, &bw[0], &v_single_mu[0]);
  InitGraph(single_mu, "Total Bandwith [kHz]", "Acceptance [%]",kRed);
  single_mu->Draw("ap");
  CMSPrelim(dataset, "#mu#tau_{h}", 0.16, 0.835);
  c1->SaveAs("thtm_single_mu.png");

  TGraph *single_mu2 = new TGraph(4, &v_single_mu_t[0], &v_single_mu[0]);
  InitGraph(single_mu2, "Threshold [GeV]", "Acceptance [%]",kRed);
  single_mu2->Draw("ap");
  CMSPrelim(dataset, "#mu#tau_{h}", 0.16, 0.835);
  c1->SaveAs("thtm_single_mu2.png");

  TGraph *mu_tau = new TGraph(4, &bw[0], &v_mu_tau[0]);
  InitGraph(mu_tau, "Total Bandwith [kHz]", "Acceptance [%]",kRed);
  mu_tau->Draw("ap");
  CMSPrelim(dataset, "#mu#tau_{h}", 0.16, 0.835);
  c1->SaveAs("thtm_mu_tau.png");

  TGraph *mu_tau2 = new TGraph(4, &v_mu_tau_t[0], &v_mu_tau[0]);
  InitGraph(mu_tau2, "#tau_{h} Threshold [GeV]", "Acceptance [%]",kRed);
  mu_tau2->Draw("ap");
  CMSPrelim(dataset, "#mu#tau_{h}", 0.16, 0.835);
  c1->SaveAs("thtm_mu_tau2.png");

  either = new TGraph(4, &bw[0], &v_either[0]);
  InitGraph(either, "Total Bandwith [kHz]", "Acceptance [%]",kBlue);
  either->Draw("ap");
  CMSPrelim(dataset, "#mu#tau_{h}", 0.16, 0.835);
  c1->SaveAs("thtm_inc.png");

  v_single_tau.clear();
  v_single_mu.clear();
  v_mu_tau.clear();
  v_either.clear();
  delete(single_tau);
  delete(single_tau2);
  delete(single_mu);
  delete(single_mu2);
  delete(mu_tau);
  delete(mu_tau2);
  delete(either);

  cout << "et ----" << endl;

  channel=float(tree->GetEntries(et));
  analysis=100*tree->GetEntries("(((t1pt>30 && abs(t1eta)<2.4)||(t2pt>30 && abs(t2eta)<2.4))&&((e1pt>20 && abs(e1eta)<2.4)||(e2pt>20 && abs(e2eta)<2.4)))*"+et)/channel;
  cout << "Analysis: " << analysis << "%" << endl;

  vector<Float_t> v_single_ele;
  vector<Float_t> v_single_ele_t;
  v_single_ele_t.push_back(31.2); v_single_ele_t.push_back(27);
  v_single_ele_t.push_back(25.1); v_single_ele_t.push_back(22.8);
  vector<Float_t> v_ele_tau;
  vector<Float_t> v_ele_tau_t;
  v_ele_tau_t.push_back(53.6); v_ele_tau_t.push_back(50);
  v_ele_tau_t.push_back(48); v_ele_tau_t.push_back(46);

  v_single_tau.push_back(100*tree->GetEntries("(triggerBits126&256)*"+et)/channel);
  v_single_tau.push_back(100*tree->GetEntries("(triggerBits180&256)*"+et)/channel);
  v_single_tau.push_back(100*tree->GetEntries("(triggerBits250&256)*"+et)/channel);
  v_single_tau.push_back(100*tree->GetEntries("(triggerBits350&256)*"+et)/channel);

  v_single_ele.push_back(100*tree->GetEntries("(triggerBits126&16)*"+et)/channel);
  v_single_ele.push_back(100*tree->GetEntries("(triggerBits180&16)*"+et)/channel);
  v_single_ele.push_back(100*tree->GetEntries("(triggerBits250&16)*"+et)/channel);
  v_single_ele.push_back(100*tree->GetEntries("(triggerBits350&16)*"+et)/channel);

  v_ele_tau.push_back(100*tree->GetEntries("(triggerBits126&1024)*"+et)/channel);
  v_ele_tau.push_back(100*tree->GetEntries("(triggerBits180&1024)*"+et)/channel);
  v_ele_tau.push_back(100*tree->GetEntries("(triggerBits250&1024)*"+et)/channel);
  v_ele_tau.push_back(100*tree->GetEntries("(triggerBits350&1024)*"+et)/channel);

  v_either.push_back(100*tree->GetEntries("((triggerBits126&256)||(triggerBits126&16)||(triggerBits126&1024))*"+et)/channel);
  v_either.push_back(100*tree->GetEntries("((triggerBits180&256)||(triggerBits180&16)||(triggerBits180&1024))*"+et)/channel);
  v_either.push_back(100*tree->GetEntries("((triggerBits250&256)||(triggerBits250&16)||(triggerBits250&1024))*"+et)/channel);
  v_either.push_back(100*tree->GetEntries("((triggerBits350&256)||(triggerBits350&16)||(triggerBits350&1024))*"+et)/channel);

  single_tau = new TGraph(4, &bw[0], &v_single_tau[0]);
  InitGraph(single_tau, "Total Bandwith [kHz]", "Acceptance [%]",kRed);
  single_tau->Draw("ap");
  CMSPrelim(dataset, "e#tau_{h}", 0.16, 0.835);
  c1->SaveAs("thte_single_tau.png");

  single_tau2 = new TGraph(4, &v_single_tau_t[0], &v_single_tau[0]);
  InitGraph(single_tau2, "Threshold [GeV]", "Acceptance [%]",kRed);
  single_tau2->Draw("ap");
  CMSPrelim(dataset, "e#tau_{h}", 0.16, 0.835);
  c1->SaveAs("thte_single_tau2.png");

  TGraph *single_ele = new TGraph(4, &bw[0], &v_single_ele[0]);
  InitGraph(single_ele, "Total Bandwith [kHz]", "Acceptance [%]",kRed);
  single_ele->Draw("ap");
  CMSPrelim(dataset, "e#tau_{h}", 0.16, 0.835);
  c1->SaveAs("thte_single_ele.png");

  TGraph *single_ele2 = new TGraph(4, &v_single_ele_t[0], &v_single_ele[0]);
  InitGraph(single_ele2, "Threshold [GeV]", "Acceptance [%]",kRed);
  single_ele2->Draw("ap");
  CMSPrelim(dataset, "e#tau_{h}", 0.16, 0.835);
  c1->SaveAs("thte_single_ele2.png");

  TGraph *ele_tau = new TGraph(4, &bw[0], &v_ele_tau[0]);
  InitGraph(ele_tau, "Total Bandwith [kHz]", "Acceptance [%]",kRed);
  ele_tau->Draw("ap");
  CMSPrelim(dataset, "e#tau_{h}", 0.16, 0.835);
  c1->SaveAs("thte_ele_tau.png");

  TGraph *ele_tau2 = new TGraph(4, &v_ele_tau_t[0], &v_ele_tau[0]);
  InitGraph(ele_tau2, "#tau_{h} Threshold [GeV]", "Acceptance [%]",kRed);
  ele_tau2->Draw("ap");
  CMSPrelim(dataset, "e#tau_{h}", 0.16, 0.835);
  c1->SaveAs("thte_ele_tau2.png");

  either = new TGraph(4, &bw[0], &v_either[0]);
  InitGraph(either, "Total Bandwith [kHz]", "Acceptance [%]",kBlue);
  either->Draw("ap");
  CMSPrelim(dataset, "e#tau_{h}", 0.16, 0.835);
  c1->SaveAs("thte_inc.png");

  v_single_tau.clear();
  v_single_ele.clear();
  v_ele_tau.clear();
  v_either.clear();
  delete(single_tau);
  delete(single_ele);
  delete(ele_tau);
  delete(single_tau2);
  delete(single_ele2);
  delete(ele_tau2);
  delete(either);

  cout << "em ----" << endl;

  channel=float(tree->GetEntries(em));
  analysis=100*tree->GetEntries("((((e1pt>20 && abs(e1eta)<2.4)||(e1pt>20 && abs(e2eta)<2.4))&&((m1pt>10 && abs(m1eta)<2.4)||(m2pt>10 && abs(m2eta)<2.4)))||(((e1pt>10 && abs(e1eta)<2.4)||(e1pt>10 && abs(e2eta)<2.4))&&((m1pt>20 && abs(m1eta)<2.4)||(m2pt>20 && abs(m2eta)<2.4))))*"+em)/channel;
  cout << "Analysis: " << analysis << "%" << endl;

  vector<Float_t> v_ele_mu;
  vector<Float_t> v_ele_mu_t;
  v_ele_mu_t.push_back(19.6); v_ele_mu_t.push_back(19.0);
  v_ele_mu_t.push_back(18.5); v_ele_mu_t.push_back(16.3);

  v_single_mu.push_back(100*tree->GetEntries("(triggerBits126&1)*"+em)/channel);
  v_single_mu.push_back(100*tree->GetEntries("(triggerBits180&1)*"+em)/channel);
  v_single_mu.push_back(100*tree->GetEntries("(triggerBits250&1)*"+em)/channel);
  v_single_mu.push_back(100*tree->GetEntries("(triggerBits350&1)*"+em)/channel);

  v_single_ele.push_back(100*tree->GetEntries("(triggerBits126&16)*"+em)/channel);
  v_single_ele.push_back(100*tree->GetEntries("(triggerBits180&16)*"+em)/channel);
  v_single_ele.push_back(100*tree->GetEntries("(triggerBits250&16)*"+em)/channel);
  v_single_ele.push_back(100*tree->GetEntries("(triggerBits350&16)*"+em)/channel);

  v_ele_mu.push_back(100*tree->GetEntries("(triggerBits126&4)*"+em)/channel);
  v_ele_mu.push_back(100*tree->GetEntries("(triggerBits180&4)*"+em)/channel);
  v_ele_mu.push_back(100*tree->GetEntries("(triggerBits250&4)*"+em)/channel);
  v_ele_mu.push_back(100*tree->GetEntries("(triggerBits350&4)*"+em)/channel);

  v_either.push_back(100*tree->GetEntries("((triggerBits126&1)||(triggerBits126&16)||(triggerBits126&4))*"+em)/channel);
  v_either.push_back(100*tree->GetEntries("((triggerBits180&1)||(triggerBits180&16)||(triggerBits180&4))*"+em)/channel);
  v_either.push_back(100*tree->GetEntries("((triggerBits250&1)||(triggerBits250&16)||(triggerBits250&4))*"+em)/channel);
  v_either.push_back(100*tree->GetEntries("((triggerBits350&1)||(triggerBits350&16)||(triggerBits350&4))*"+em)/channel);

  single_mu = new TGraph(4, &bw[0], &v_single_mu[0]);
  InitGraph(single_mu, "Total Bandwith [kHz]", "Acceptance [%]",kRed);
  single_mu->Draw("ap");
  CMSPrelim(dataset, "e#mu", 0.16, 0.835);
  c1->SaveAs("tmte_single_mu.png");

  single_mu2 = new TGraph(4, &v_single_mu_t[0], &v_single_mu[0]);
  InitGraph(single_mu2, "Threshold [GeV]", "Acceptance [%]",kRed);
  single_mu2->Draw("ap");
  CMSPrelim(dataset, "e#mu", 0.16, 0.835);
  c1->SaveAs("tmte_single_mu2.png");

  single_ele = new TGraph(4, &bw[0], &v_single_ele[0]);
  InitGraph(single_ele, "Total Bandwith [kHz]", "Acceptance [%]",kRed);
  single_ele->Draw("ap");
  CMSPrelim(dataset, "e#mu", 0.16, 0.835);
  c1->SaveAs("tmte_single_ele.png");

  single_ele2 = new TGraph(4, &v_single_ele_t[0], &v_single_ele[0]);
  InitGraph(single_ele2, "Threshold [GeV]", "Acceptance [%]",kRed);
  single_ele2->Draw("ap");
  CMSPrelim(dataset, "e#mu", 0.16, 0.835);
  c1->SaveAs("tmte_single_ele2.png");

  TGraph *ele_mu = new TGraph(4, &bw[0], &v_ele_mu[0]);
  InitGraph(ele_mu, "Total Bandwith [kHz]", "Acceptance [%]",kRed);
  ele_mu->Draw("ap");
  CMSPrelim(dataset, "e#mu", 0.16, 0.835);
  c1->SaveAs("tmte_ele_mu.png");

  TGraph *ele_mu2 = new TGraph(4, &v_ele_mu_t[0], &v_ele_mu[0]);
  InitGraph(ele_mu2, "e Threshold [GeV]", "Acceptance [%]",kRed);
  ele_mu2->Draw("ap");
  CMSPrelim(dataset, "e#mu", 0.16, 0.835);
  c1->SaveAs("tmte_ele_mu2.png");

  either = new TGraph(4, &bw[0], &v_either[0]);
  InitGraph(either, "Total Bandwith [kHz]", "Acceptance [%]",kBlue);
  either->Draw("ap");
  CMSPrelim(dataset, "e#mu", 0.16, 0.835);
  c1->SaveAs("tmte_inc.png");

}
