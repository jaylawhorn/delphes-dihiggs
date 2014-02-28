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

void test_iap(const TString inputfile="root://eoscms.cern.ch//store/group/phys_higgs/upgrade/PhaseII/Configuration4v2/140PileUp/ttB-4p-0-900-v1510_14TEV/ttB-4p-0-900-v1510_14TEV_100133601_PhaseII_Conf4v2_140PileUp.root", 
		   const Float_t xsec=2.6673,
		   const TString outputfile="test.root") {

  // read input input file
  TChain chain("Delphes");
  chain.Add(inputfile);
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  TClonesArray *branchEvent = treeReader->UseBranch("Event");

  LHEFEvent *event;

  // set up output variables and file
  Int_t nEvents;
  Float_t eventWeight;

  TFile *outFile = new TFile(outputfile, "RECREATE");

  // tree to hold the number of events in the file before selection
  TTree *sampTree = new TTree("Info", "Info");
  sampTree->Branch("nEvents",       &nEvents,        "nEvents/i");
  nEvents=numberOfEntries;
  sampTree->Fill();

  // tree to hold information about selected events
  TTree *outTree = new TTree("Events", "Events");
  outTree->Branch("eventWeight",    &eventWeight,    "eventWeight/f");  // event weight from cross-section and Event->Weight

  for (Int_t iEntry=0; iEntry<numberOfEntries; iEntry++) { // entry loop
    treeReader->ReadEntry(iEntry);

    // comment out following line for di-higgs samples
    event = (LHEFEvent*) branchEvent->At(0);
    eventWeight = 1;
    eventWeight *= xsec;
    // comment out following line for di-higgs samples
    eventWeight *= event->Weight;

    outTree->Fill();

  } // end event loop

  outFile->Write();
  outFile->Save();

}
