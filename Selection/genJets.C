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

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;

Float_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 );

void genJets(const TString inputfile="root://eoscms.cern.ch//store/group/phys_higgs/upgrade/PhaseII/Configuration4v2/140PileUp/HHToTTBB_14TeV/HHToTTBB_14TeV_1.root", 
		   const Float_t xsec=2.92,
		   const TString outputfile="test.root") {

  const Int_t B_ID_CODE = 5;
  const Int_t TAU_ID_CODE = 15;

  // read input input file
  TChain chain("Delphes");
  chain.Add(inputfile);
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  TClonesArray *branchParticle = treeReader->UseBranch("Particle");

  //set up loop variables
  GenParticle *genParticle;

  // set up storage variables
  GenParticle *gentau1, *gentau2;
  GenParticle *genb1, *genb2;

  Int_t iGenTau1=-1,    iGenTau2=-1;
  Int_t iGenB1=-1,      iGenB2=-1;

  // set up output variables and file
  Int_t nEvents;
  Float_t eventWeight;
  LorentzVector *genTau1=0, *genTau2=0;
  LorentzVector *genB1=0, *genB2=0;

  TFile *outFile = new TFile(outputfile, "RECREATE");

  // tree to hold the number of events in the file before selection
  TTree *sampTree = new TTree("Info", "Info");
  sampTree->Branch("nEvents",       &nEvents,        "nEvents/i");
  nEvents=numberOfEntries;
  sampTree->Fill();

  cout << nEvents << endl;

  // tree to hold information about selected events
  TTree *outTree = new TTree("Events", "Events");
  outTree->Branch("eventWeight",    &eventWeight,    "eventWeight/f");  // event weight from cross-section and Event->Weight
  outTree->Branch("genTau1",        "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &genTau1);      // 4-vector for reconstructed leading tau
  outTree->Branch("genTau2",        "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &genTau2);      // 4-vector for reconstructed second tau
  outTree->Branch("genB1",          "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &genB1);        // 4-vector for reconstructed leading b-jet
  outTree->Branch("genB2",          "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &genB2);        // 4-vector for reconstructed second b-jet

  // define placeholder vector for things that don't exist
  LorentzVector nothing(999,999,0,999);

  for (Int_t iEntry=0; iEntry<numberOfEntries; iEntry++) { // entry loop
    treeReader->ReadEntry(iEntry);

    iGenTau1=-1;    iGenTau2=-1;
    iGenB1=-1;      iGenB2=-1;

    eventWeight = 1;
    eventWeight *= xsec;

    for (Int_t iParticle=0; iParticle<branchParticle->GetEntries(); iParticle++) { // generator particle loop
      genParticle = (GenParticle*) branchParticle->At(iParticle);

      if ( fabs(genParticle->PID) == TAU_ID_CODE ) { // tau switch
        if (iGenTau1==-1) {
          iGenTau1 = iParticle;
          gentau1 = (GenParticle*) branchParticle->At(iGenTau1);
        }
        else if ((iGenTau1!=-1) && ( genParticle->PT > gentau1->PT)) {
          iGenTau2 = iGenTau1;
          gentau2 = (GenParticle*) branchParticle->At(iGenTau2);
          iGenTau1 = iParticle;
          gentau1 = (GenParticle*) branchParticle->At(iGenTau1);
        }
	else if ((iGenTau1!=-1) && (iGenTau2==-1)) {
          iGenTau2 = iParticle;
          gentau2 = (GenParticle*) branchParticle->At(iGenTau2);
        }
	else if ((iGenTau1!=-1) && (genParticle->PT > gentau2->PT)) {
	  iGenTau2 = iParticle;
          gentau2 = (GenParticle*) branchParticle->At(iGenTau2);
	}
      }

      if ( fabs(genParticle->PID) == B_ID_CODE ) { // b-quark switch
        if (iGenB1==-1) {
          iGenB1 = iParticle;
          genb1 = (GenParticle*) branchParticle->At(iGenB1);
        }
        else if ((iGenB1!=-1) && ( genParticle->PT > genb1->PT)) {
          iGenB2 = iGenB1;
          genb2 = (GenParticle*) branchParticle->At(iGenB2);
          iGenB1 = iParticle;
          genb1 = (GenParticle*) branchParticle->At(iGenB1);
        }
	else if ((iGenB1!=-1) && (iGenB2==-1)) {
          iGenB2 = iParticle;
          genb2 = (GenParticle*) branchParticle->At(iGenB2);
        }
	else if ((iGenB1!=-1) && (genParticle->PT > genb2->PT)) {
	  iGenB2 = iParticle;
          genb2 = (GenParticle*) branchParticle->At(iGenB2);
	}
      }

    } // end particle loop

    cout << iGenB1 << " " << iGenB2 << endl;

    // store generator particles
    LorentzVector vGenTau1(0,0,0,0);
    if ( iGenTau1 == -1) genTau1 = &nothing;
    else {
      vGenTau1.SetPt(gentau1->PT);
      vGenTau1.SetEta(gentau1->Eta);
      vGenTau1.SetPhi(gentau1->Phi);
      vGenTau1.SetM(gentau1->Mass);
      genTau1 = &vGenTau1;
    }

    LorentzVector vGenTau2(0,0,0,0);
    if ( iGenTau2 == -1) genTau2 = &nothing;
    else {
      vGenTau2.SetPt(gentau2->PT);
      vGenTau2.SetEta(gentau2->Eta);
      vGenTau2.SetPhi(gentau2->Phi);
      vGenTau2.SetM(gentau2->Mass);
      genTau2 = &vGenTau2;
    }

    LorentzVector vGenB1(0,0,0,0);
    if ( iGenB1 == -1) genB1 = &nothing;
    else {
      vGenB1.SetPt(genb1->PT);
      vGenB1.SetEta(genb1->Eta);
      vGenB1.SetPhi(genb1->Phi);
      vGenB1.SetM(genb1->Mass);
      genB1 = &vGenB1;
    }

    LorentzVector vGenB2(0,0,0,0);
    if ( iGenB2 == -1) genB2 = &nothing;
    else {
      vGenB2.SetPt(genb2->PT);
      vGenB2.SetEta(genb2->Eta);
      vGenB2.SetPhi(genb2->Phi);
      vGenB2.SetM(genb2->Mass);
      genB2 = &vGenB2;
    }

    outTree->Fill();

  } // end event loop

  outFile->Write();
  outFile->Save();

}

Float_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 ) {

  const Float_t pi = 3.14159265358979;

  Float_t etaDiff = (eta1-eta2);
  Float_t phiDiff = fabs(phi1-phi2);
  while ( phiDiff>pi ) phiDiff = fabs(phiDiff-2.0*pi);

  Float_t deltaRSquared = etaDiff*etaDiff + phiDiff*phiDiff;

  return TMath::Sqrt(deltaRSquared);

}
