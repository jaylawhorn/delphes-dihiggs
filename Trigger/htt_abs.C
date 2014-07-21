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

void htt_abs() {

  TFile *file = new TFile("/afs/cern.ch/work/j/jlawhorn/htt.root");
  TTree *tree = (TTree*)file->Get("Events");
  Trigger::Event data;
  tree->SetBranchAddress("Events",&data);

  Float_t tot=tree->GetEntries();

  TString tt="t1pt>0 && t2pt>0";
  TString mt="m1pt>0 && (t1pt>0 || t2pt>0)";
  TString et="e1pt>0 && (t1pt>0 || t2pt>0)";
  TString em="m1pt>0 && e1pt>0";
  TString ee="e1pt>0 && e2pt>0";
  TString mm="m1pt>0 && m2pt>0";
  
  cout << "Total entries: " << tot << endl;
  cout << "Tau-tau:       " << 100*tree->GetEntries(tt)/tot << "%" << endl;
  cout << "Mu-tau:        " << 100*tree->GetEntries(mt)/tot << "%" << endl;
  cout << "Ele-tau:       " << 100*tree->GetEntries(et)/tot << "%" << endl;
  cout << "Ele-mu:        " << 100*tree->GetEntries(em)/tot << "%" << endl;
  cout << "ee:            " << 100*tree->GetEntries(ee)/tot << "%" << endl;
  cout << "mm:            " << 100*tree->GetEntries(mm)/tot << "%" << endl;
  cout << "Nothing:        " << 100*tree->GetEntries("!("+tt+")*!("+mt+")*!("+et+")*!("+em+")*!("+ee+")*!("+mm+")")/tot << "%" << endl;
  cout << endl;

  cout << "tt ----" << endl;
  cout << "      analysis: " << 100*tree->GetEntries("(t1pt>45 && t2pt>45 && abs(t1eta)<2.4 && abs(t2eta)<2.4)*"+tt)/tree->GetEntries(tt) << "%" << endl;
  cout << "    single tau: " << 100*tree->GetEntries("(triggerBits180&256)*"+tt)/tree->GetEntries(tt) << "%" << endl;
  cout << "    double tau: " << 100*tree->GetEntries("(triggerBits180&512)*"+tt)/tree->GetEntries(tt) << "%" << endl;
  cout << "          both: " << 100*tree->GetEntries("(triggerBits180&512)*(triggerBits180&256)*"+tt)/tree->GetEntries(tt) << "%" << endl;

  cout << "mt ----" << endl;
  cout << "      analysis: " << 100*tree->GetEntries("(t1pt>30 && m1pt>20 && abs(t1eta)<2.4 && abs(m1eta)<2.4)*"+mt)/tree->GetEntries(mt) << "%" << endl;
  cout << "    single tau: " << 100*tree->GetEntries("(triggerBits180&256)*"+mt)/tree->GetEntries(mt) << "%" << endl;
  cout << "    single mu : " << 100*tree->GetEntries("(triggerBits180&1)*"+mt)/tree->GetEntries(mt) << "%" << endl;
  cout << "    mu-tau    : " << 100*tree->GetEntries("(triggerBits180&2048)*"+mt)/tree->GetEntries(mt) << "%" << endl;
  cout << "      tau + mu: " << 100*tree->GetEntries("(triggerBits180&256)*(triggerBits180&1)*"+mt)/tree->GetEntries(mt) << "%" << endl;
  cout << "  tau + mu-tau: " << 100*tree->GetEntries("(triggerBits180&256)*(triggerBits180&2048)*"+mt)/tree->GetEntries(mt) << "%" << endl;
  cout << "   mu + mu-tau: " << 100*tree->GetEntries("(triggerBits180&1)*(triggerBits180&2048)*"+mt)/tree->GetEntries(mt) << "%" << endl;

  cout << "et ----" << endl;
  cout << "      analysis: " << 100*tree->GetEntries("(t1pt>30 && e1pt>20 && abs(t1eta)<2.4 && abs(e1eta)<2.4)*"+et)/tree->GetEntries(et) << "%" << endl;
  cout << "    single tau: " << 100*tree->GetEntries("(triggerBits180&256)*"+et)/tree->GetEntries(et) << "%" << endl;
  cout << "    single ele: " << 100*tree->GetEntries("(triggerBits180&16)*"+et)/tree->GetEntries(et) << "%" << endl;
  cout << "    ele-tau   : " << 100*tree->GetEntries("(triggerBits180&1024)*"+et)/tree->GetEntries(et) << "%" << endl;
  cout << "     tau + ele: " << 100*tree->GetEntries("(triggerBits180&256)*(triggerBits180&16)*"+et)/tree->GetEntries(et) << "%" << endl;
  cout << " tau + ele-tau: " << 100*tree->GetEntries("(triggerBits180&256)*(triggerBits180&1024)*"+et)/tree->GetEntries(et) << "%" << endl;
  cout << " ele + ele-tau: " << 100*tree->GetEntries("(triggerBits180&16)*(triggerBits180&1024)*"+et)/tree->GetEntries(et) << "%" << endl;

  cout << "em ----" << endl;
  cout << "      analysis: " << 100*tree->GetEntries("(e1pt>20 && m1pt>20 && abs(e1eta)<2.4 && abs(e1eta)<2.4)*"+em)/tree->GetEntries(em) << "%" << endl;
  cout << "    single ele: " << 100*tree->GetEntries("(triggerBits180&8)*"+em)/tree->GetEntries(em) << "%" << endl;
  cout << "    single mu : " << 100*tree->GetEntries("(triggerBits180&1)*"+em)/tree->GetEntries(em) << "%" << endl;
  cout << "    ele-mu    : " << 100*tree->GetEntries("(triggerBits180&4)*"+em)/tree->GetEntries(em) << "%" << endl;
  cout << "      ele + mu: " << 100*tree->GetEntries("(triggerBits180&8)*(triggerBits180&1)*"+em)/tree->GetEntries(em) << "%" << endl;
  cout << "  ele + ele-mu: " << 100*tree->GetEntries("(triggerBits180&8)*(triggerBits180&4)*"+em)/tree->GetEntries(em) << "%" << endl;
  cout << "   mu + ele-mu: " << 100*tree->GetEntries("(triggerBits180&4)*(triggerBits180&1)*"+em)/tree->GetEntries(em) << "%" << endl;


  /*
  cout << " ----" << endl;
  cout << "tt/mt " << 100*tree->GetEntries("("+tt+")*("+mt+")")/tot << endl;
  cout << "tt/et " << 100*tree->GetEntries("("+tt+")*("+et+")")/tot << endl;
  cout << "tt/em " << 100*tree->GetEntries("("+tt+")*("+em+")")/tot << endl;
  cout << "tt/ee " << 100*tree->GetEntries("("+tt+")*("+ee+")")/tot << endl;
  cout << "tt/mm " << 100*tree->GetEntries("("+tt+")*("+mm+")")/tot << endl;
  cout << " ----" << endl;
  cout << "mt/et " << 100*tree->GetEntries("("+mt+")*("+et+")")/tot << endl;
  cout << "mt/em " << 100*tree->GetEntries("("+mt+")*("+em+")")/tot << endl;
  cout << "mt/ee " << 100*tree->GetEntries("("+mt+")*("+ee+")")/tot << endl;
  cout << "mt/mm " << 100*tree->GetEntries("("+mt+")*("+mm+")")/tot << endl;
  cout << " ----" << endl;
  cout << "et/em " << 100*tree->GetEntries("("+et+")*("+em+")")/tot << endl;
  cout << "et/ee " << 100*tree->GetEntries("("+et+")*("+ee+")")/tot << endl;
  cout << "et/mm " << 100*tree->GetEntries("("+et+")*("+mm+")")/tot << endl;
  cout << " ----" << endl;
  cout << "em/ee " << 100*tree->GetEntries("("+em+")*("+ee+")")/tot << endl;
  cout << "em/mm " << 100*tree->GetEntries("("+em+")*("+mm+")")/tot << endl;
  cout << " ----" << endl;
  cout << "ee/mm " << 100*tree->GetEntries("("+ee+")*("+mm+")")/tot << endl;

  cout << " ----" << endl;
  cout << "tt/mt/et " << 100*tree->GetEntries("("+tt+")*("+mt+")*("+et+")")/tot << endl;
  cout << "tt/mt/em " << 100*tree->GetEntries("("+tt+")*("+mt+")*("+em+")")/tot << endl;
  cout << "tt/mt/ee " << 100*tree->GetEntries("("+tt+")*("+mt+")*("+ee+")")/tot << endl;
  cout << "tt/mt/mm " << 100*tree->GetEntries("("+tt+")*("+mt+")*("+mm+")")/tot << endl;
  cout << "tt/et/em " << 100*tree->GetEntries("("+tt+")*("+et+")*("+em+")")/tot << endl;
  cout << "tt/et/ee " << 100*tree->GetEntries("("+tt+")*("+et+")*("+ee+")")/tot << endl;
  cout << "tt/et/mm " << 100*tree->GetEntries("("+tt+")*("+et+")*("+mm+")")/tot << endl;
  cout << "tt/em/ee " << 100*tree->GetEntries("("+tt+")*("+em+")*("+ee+")")/tot << endl;
  cout << "tt/em/mm " << 100*tree->GetEntries("("+tt+")*("+em+")*("+mm+")")/tot << endl;
  cout << "tt/ee/mm " << 100*tree->GetEntries("("+tt+")*("+ee+")*("+mm+")")/tot << endl;
  cout << " ----" << endl;
  cout << "mt/et/em " << 100*tree->GetEntries("("+mt+")*("+et+")*("+em+")")/tot << endl;
  cout << "mt/et/ee " << 100*tree->GetEntries("("+mt+")*("+et+")*("+ee+")")/tot << endl;
  cout << "mt/et/mm " << 100*tree->GetEntries("("+mt+")*("+et+")*("+mm+")")/tot << endl;
  cout << "mt/em/ee " << 100*tree->GetEntries("("+mt+")*("+em+")*("+ee+")")/tot << endl;
  cout << "mt/em/mm " << 100*tree->GetEntries("("+mt+")*("+em+")*("+mm+")")/tot << endl;
  cout << "mt/ee/mm " << 100*tree->GetEntries("("+mt+")*("+ee+")*("+mm+")")/tot << endl;
  cout << " ----" << endl;
  cout << "et/em/ee " << 100*tree->GetEntries("("+et+")*("+em+")*("+ee+")")/tot << endl;
  cout << "et/em/mm " << 100*tree->GetEntries("("+et+")*("+em+")*("+mm+")")/tot << endl;
  cout << " ----" << endl;
  cout << "em/ee/mm " << 100*tree->GetEntries("("+em+")*("+ee+")*("+mm+")")/tot << endl;

  cout << " ----" << endl;
  cout << "tt/mt/et/em " << 100*tree->GetEntries("("+tt+")*("+mt+")*("+et+")*("+em+")")/tot << endl;
  cout << "tt/mt/et/ee " << 100*tree->GetEntries("("+tt+")*("+mt+")*("+et+")*("+ee+")")/tot << endl;
  cout << "tt/mt/et/mm " << 100*tree->GetEntries("("+tt+")*("+mt+")*("+et+")*("+mm+")")/tot << endl;
  cout << "tt/mt/ee/mm " << 100*tree->GetEntries("("+tt+")*("+mt+")*("+ee+")*("+mm+")")/tot << endl;
  cout << "tt/et/em/ee " << 100*tree->GetEntries("("+tt+")*("+et+")*("+em+")*("+ee+")")/tot << endl;
  cout << "tt/et/em/mm " << 100*tree->GetEntries("("+tt+")*("+et+")*("+em+")*("+mm+")")/tot << endl;
  cout << "tt/et/ee/mm " << 100*tree->GetEntries("("+tt+")*("+et+")*("+ee+")*("+mm+")")/tot << endl;
  cout << "tt/em/ee/mm " << 100*tree->GetEntries("("+tt+")*("+em+")*("+ee+")*("+mm+")")/tot << endl;
  cout << " ----" << endl;
  cout << "mt/et/em/ee " << 100*tree->GetEntries("("+mt+")*("+et+")*("+em+")*("+ee+")")/tot << endl;
  cout << "mt/et/em/mm " << 100*tree->GetEntries("("+mt+")*("+et+")*("+em+")*("+mm+")")/tot << endl;
  cout << "mt/et/ee/mm " << 100*tree->GetEntries("("+mt+")*("+et+")*("+ee+")*("+mm+")")/tot << endl;
  cout << "mt/em/ee/mm " << 100*tree->GetEntries("("+mt+")*("+em+")*("+ee+")*("+mm+")")/tot << endl;
  cout << " ----" << endl;
  cout << "et/em/ee/mm " << 100*tree->GetEntries("("+et+")*("+em+")*("+ee+")*("+mm+")")/tot << endl;
  cout << " ----" << endl;
  cout << "tt/mt/et/em/ee/mm " << 100*tree->GetEntries("("+tt+")*("+mt+")*("+et+")*("+em+")*("+ee+")*("+mm+")")/tot << endl;
  */
}
