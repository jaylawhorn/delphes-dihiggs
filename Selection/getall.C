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
#include "mt2.hh"
#endif

void getall(const Float_t xsec=2.92) { 

  char inputfile[300];

  TChain chain("Delphes");

  for (Int_t i=0; i<175; i++) {
    sprintf(inputfile,"root://eoscms.cern.ch//store/group/phys_higgs/upgrade/PhaseII/Configuration4v2/140PileUp/HHToTTBB_14TeV/HHToTTBB_14TeV_%i.root",i);
    chain.Add(inputfile);
  }
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  Int_t nEvents;
  Double_t eventWeight;
  Double_t totalWeight;

  for (Int_t iEntry=0; iEntry<numberOfEntries; iEntry++) { // entry loop
    treeReader->ReadEntry(iEntry);

    // comment out following line for di-higgs samples
    //event = (LHEFEvent*) branchEvent->At(0);
    eventWeight = 1;
    eventWeight *= xsec;
    // comment out following line for di-higgs samples
    //eventWeight *= event->Weight;
    totalWeight+=eventWeight;
    nEvents+=1;
  }
  cout << numberOfEntries << " " << nEvents << endl;
  cout << "Total cross section: " << totalWeight/float(nEvents) << endl;
  cout << "At 3000/fb         : " << 3000*totalWeight/float(nEvents) << endl;

}
