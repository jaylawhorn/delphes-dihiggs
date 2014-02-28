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

//#include "bJetScaleCorr.h"
#endif

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;

Float_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 );

void testWeighting(const TString inputfile="root://eoscms.cern.ch//store/group/phys_higgs/upgrade/PhaseII/Configuration3/140PileUp/BB-4p-0-300-v1510_14TEV/BB-4p-0-300-v1510_14TEV_99808731_PhaseII_Conf3_140PileUp.root",
		   const Float_t xsec=249.97710,
		   const TString outputfile="test.root") {

  // read input input file
  TChain chain("Delphes");
  chain.Add(inputfile);
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  TClonesArray *branchEvent = treeReader->UseBranch("Event");

  LHEFEvent *event;

  Int_t nEvents;
  Float_t eventWeight;

  TFile *outFile = new TFile(outputfile, "RECREATE");

  TTree *sampTree = new TTree("Info", "Info");
  sampTree->Branch("nEvents",       &nEvents,        "nEvents/i");
  nEvents=numberOfEntries;
  sampTree->Fill();

  TTree *outTree = new TTree("Events", "Events");
  outTree->Branch("eventWeight",    &eventWeight,    "eventWeight/f");

  for (Int_t iEntry=0; iEntry<numberOfEntries; iEntry++) { // entry loop
    treeReader->ReadEntry(iEntry);

    //cout << " ---- " << endl;

    event = (LHEFEvent*) branchEvent->At(0);
    //cout << xsec << ", " << numberOfEntries << ", " << event->Weight << ", ";
    eventWeight = 1;
    eventWeight = xsec;
    eventWeight *= event->Weight;
    //cout << eventWeight << endl;

    outTree->Fill();

  } // end event loop

  outFile->Write();
  outFile->Close();

  cout << "----SUMMARY----" << endl;
  cout << " input file " << inputfile << " selection done " << endl;

}

Float_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 ) {

  const Float_t pi = 3.14159265358979;

  Float_t etaDiff = (eta1-eta2);
  Float_t phiDiff = fabs(phi1-phi2);
  while ( phiDiff>pi ) phiDiff = fabs(phiDiff-2.0*pi);

  Float_t deltaRSquared = etaDiff*etaDiff + phiDiff*phiDiff;

  return TMath::Sqrt(deltaRSquared);

}
