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
#include "Math/LorentzVector.h"
#include "modules/Delphes.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "classes/DelphesClasses.h"

#endif

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;

void testNewSelect(const TString infile1="/afs/cern.ch/work/k/klawhorn/SnowmassSamples/PhaseII/Configuration4v2/HHToTTBB_14TeV_1.root",
		   const TString infile2="test.root") {

  // set up output variables and file
  LorentzVector *genTau1=0, *genTau2=0, *genDecayTau1=0, *genDecayTau2=0, *recoTau1=0, *recoTau2=0;
  LorentzVector *genB1=0, *genB2=0, *recoB1=0, *recoB2=0;

  TFile *file1= new TFile(infile1, "READ"); assert(file1);
  TFile *file2= new TFile(infile2, "READ"); assert(file2);

  TTree *tree = (TTree*) file1->Get("Events"); assert(tree);
  tree->SetBranchAddress("genTau1",        &genTau1);
  tree->SetBranchAddress("genTau2",        &genTau2);
  tree->SetBranchAddress("genDecayTau1",   &genDecayTau1);
  tree->SetBranchAddress("genDecayTau2",   &genDecayTau2);
  tree->SetBranchAddress("recoTau1",       &recoTau1);
  tree->SetBranchAddress("recoTau2",       &recoTau2);
  tree->SetBranchAddress("genB1",          &genB1);
  tree->SetBranchAddress("genB2",          &genB2);
  tree->SetBranchAddress("recoB1",         &recoB1);
  tree->SetBranchAddress("recoB2",         &recoB2);

  TH1D *hold = new TH1D("old", "old", 100, 0, 300);
  TH1D *hnew = new TH1D("new", "new", 100, 0, 300);

  for (Int_t iEntry=0; iEntry<tree->GetEntries(); iEntry++) { // entry loop
    //cout << recoTau1->Pt() << endl;
    tree->GetEntry(iEntry);
    hold->Fill(recoB1->Pt());
    hold->Fill(recoB2->Pt());
  } // end event loop

  tree = (TTree*) file2->Get("Events"); assert(tree);
  tree->SetBranchAddress("genTau1",        &genTau1);
  tree->SetBranchAddress("genTau2",        &genTau2);
  tree->SetBranchAddress("genDecayTau1",   &genDecayTau1);
  tree->SetBranchAddress("genDecayTau2",   &genDecayTau2);
  tree->SetBranchAddress("recoTau1",       &recoTau1);
  tree->SetBranchAddress("recoTau2",       &recoTau2);
  tree->SetBranchAddress("genB1",          &genB1);
  tree->SetBranchAddress("genB2",          &genB2);
  tree->SetBranchAddress("recoB1",         &recoB1);
  tree->SetBranchAddress("recoB2",         &recoB2);

  for (Int_t iEntry=0; iEntry<tree->GetEntries(); iEntry++) { // entry loop
    tree->GetEntry(iEntry);
    hnew->Fill(recoB1->Pt());
    hnew->Fill(recoB2->Pt());
  } // end event loop

  TCanvas *c = new TCanvas("c", "c", 800, 600);
  hold->Draw();
  hnew->SetLineColor(kRed);
  hnew->Draw("same");

}
