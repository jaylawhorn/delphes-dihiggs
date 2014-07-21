#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TVector2.h>
#include <TMath.h>
#include <TChain.h>
#include <TH1.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include "Math/LorentzVector.h"
#include "TLorentzVector.h"

#include "triggerMenu.h"

#endif

void whbb_abs() {

  TFile *file = new TFile("/afs/cern.ch/work/j/jlawhorn/whbb.root");
  TTree *tree = (TTree*)file->Get("Events");
  Trigger::Event data;
  tree->SetBranchAddress("Events",&data);

  TString we="(e1pt>0 || e2pt>0)";
  TString wm="(m1pt>0 || m2pt>0)";
  TString wt="(t1pt>0 || t2pt>0)";

  Float_t tot=tree->GetEntries();

  cout << "Total entries: " << tot << endl;

  cout << "W->enu:    " << 100*tree->GetEntries(we)/tot << "%" << endl;
  cout << "W->munu:   " << 100*tree->GetEntries(wm)/tot << "%" << endl;
  cout << "W->taunu:  " << 100*tree->GetEntries(wt)/tot << "%" << endl;
  cout << "(me)   " << 100*tree->GetEntries(wm+"*"+we)/tot << "%" << endl;
  cout << "(mt)   " << 100*tree->GetEntries(wm+"*"+wt)/tot << "%" << endl;
  cout << "(et)   " << 100*tree->GetEntries(we+"*"+wt)/tot << "%" << endl;
  cout << "(emt)  " << 100*tree->GetEntries(we+"*"+wt+"*"+wm)/tot << "%" << endl;
  cout << endl;

  cout << "** 126kHz menu" << endl;
  cout << endl;

  cout << "W->enu ---" << endl;
  cout <<     "kSingleEle: " << 100*tree->GetEntries("(triggerBits126&16)*"+we)/tree->GetEntries(we) << "%" << endl;
  cout <<     "kDoubleJet: " << 100*tree->GetEntries("(triggerBits126&8192)*"+we)/tree->GetEntries(we) << "%" << endl;
  cout <<     "kEleJet   : " << 100*tree->GetEntries("(triggerBits126&32768)*"+we)/tree->GetEntries(we) << "%" << endl;
  cout <<     "kEleMet   : " << 100*tree->GetEntries("(triggerBits126&131072)*"+we)/tree->GetEntries(we) << "%" << endl;
  cout << "W->munu --" << endl;
  cout <<     "kSingleMu : " << 100*tree->GetEntries("(triggerBits126&1)*"+wm)/tree->GetEntries(wm) << "%" << endl;
  cout <<     "kDoubleJet: " << 100*tree->GetEntries("(triggerBits126&8192)*"+wm)/tree->GetEntries(wm) << "%" << endl;
  cout <<     "kMuJet    : " << 100*tree->GetEntries("(triggerBits126&65536)*"+wm)/tree->GetEntries(wm) << "%" << endl;
  cout <<     "kMuMet    : " << 100*tree->GetEntries("(triggerBits126&262144)*"+wm)/tree->GetEntries(wm) << "%" << endl;

  cout << endl;
  cout << "** 180kHz menu" << endl;
  cout << endl;

  cout << "W->enu ---" << endl;
  cout <<     "kSingleEle: " << 100*tree->GetEntries("(triggerBits180&16)*"+we)/tree->GetEntries(we) << "%" << endl;
  cout <<     "kDoubleJet: " << 100*tree->GetEntries("(triggerBits180&8192)*"+we)/tree->GetEntries(we) << "%" << endl;
  cout <<     "kEleJet   : " << 100*tree->GetEntries("(triggerBits180&32768)*"+we)/tree->GetEntries(we) << "%" << endl;
  cout <<     "kEleMet   : " << 100*tree->GetEntries("(triggerBits180&131072)*"+we)/tree->GetEntries(we) << "%" << endl;
  cout << "W->munu --" << endl;
  cout <<     "kSingleMu : " << 100*tree->GetEntries("(triggerBits180&1)*"+wm)/tree->GetEntries(wm) << "%" << endl;
  cout <<     "kDoubleJet: " << 100*tree->GetEntries("(triggerBits180&8192)*"+wm)/tree->GetEntries(wm) << "%" << endl;
  cout <<     "kMuJet    : " << 100*tree->GetEntries("(triggerBits180&65536)*"+wm)/tree->GetEntries(wm) << "%" << endl;
  cout <<     "kMuMet    : " << 100*tree->GetEntries("(triggerBits180&262144)*"+wm)/tree->GetEntries(wm) << "%" << endl;

  cout << endl;
  cout << "** 250kHz menu" << endl;
  cout << endl;

  cout << "W->enu ---" << endl;
  cout <<     "kSingleEle: " << 100*tree->GetEntries("(triggerBits250&16)*"+we)/tree->GetEntries(we) << "%" << endl;
  cout <<     "kDoubleJet: " << 100*tree->GetEntries("(triggerBits250&8192)*"+we)/tree->GetEntries(we) << "%" << endl;
  cout <<     "kEleJet   : " << 100*tree->GetEntries("(triggerBits250&32768)*"+we)/tree->GetEntries(we) << "%" << endl;
  cout <<     "kEleMet   : " << 100*tree->GetEntries("(triggerBits250&131072)*"+we)/tree->GetEntries(we) << "%" << endl;
  cout << "W->munu --" << endl;
  cout <<     "kSingleMu : " << 100*tree->GetEntries("(triggerBits250&1)*"+wm)/tree->GetEntries(wm) << "%" << endl;
  cout <<     "kDoubleJet: " << 100*tree->GetEntries("(triggerBits250&8192)*"+wm)/tree->GetEntries(wm) << "%" << endl;
  cout <<     "kMuJet    : " << 100*tree->GetEntries("(triggerBits250&65536)*"+wm)/tree->GetEntries(wm) << "%" << endl;
  cout <<     "kMuMet    : " << 100*tree->GetEntries("(triggerBits250&262144)*"+wm)/tree->GetEntries(wm) << "%" << endl;

  cout << endl;
  cout << "** 350kHz menu" << endl;
  cout << endl;

  cout << "W->enu ---" << endl;
  cout <<     "kSingleEle: " << 100*tree->GetEntries("(triggerBits350&16)*"+we)/tree->GetEntries(we) << "%" << endl;
  cout <<     "kDoubleJet: " << 100*tree->GetEntries("(triggerBits350&8192)*"+we)/tree->GetEntries(we) << "%" << endl;
  cout <<     "kEleJet   : " << 100*tree->GetEntries("(triggerBits350&32768)*"+we)/tree->GetEntries(we) << "%" << endl;
  cout <<     "kEleMet   : " << 100*tree->GetEntries("(triggerBits350&131072)*"+we)/tree->GetEntries(we) << "%" << endl;
  cout << "W->munu --" << endl;
  cout <<     "kSingleMu : " << 100*tree->GetEntries("(triggerBits350&1)*"+wm)/tree->GetEntries(wm) << "%" << endl;
  cout <<     "kDoubleJet: " << 100*tree->GetEntries("(triggerBits350&8192)*"+wm)/tree->GetEntries(wm) << "%" << endl;
  cout <<     "kMuJet    : " << 100*tree->GetEntries("(triggerBits350&65536)*"+wm)/tree->GetEntries(wm) << "%" << endl;
  cout <<     "kMuMet    : " << 100*tree->GetEntries("(triggerBits350&262144)*"+wm)/tree->GetEntries(wm) << "%" << endl;

}
