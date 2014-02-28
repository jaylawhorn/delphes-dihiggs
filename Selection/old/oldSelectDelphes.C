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

void selectDelphes(const TString inputfile="root://eoscms.cern.ch//store/group/phys_higgs/upgrade/PhaseII/Configuration4v2/140PileUp/HHToTTBB_14TeV/HHToTTBB_14TeV_1.root", 
		   const Float_t xsec=2.92,
		   const TString outputfile="test2.root") {

  // declare constants
  const Double_t MUON_MASS = 0.105658369;
  const Double_t ELE_MASS  = 0.000511;

  const Int_t B_ID_CODE = 5;
  const Int_t TAU_ID_CODE = 15;

  const Float_t MAX_MATCH_DIST = 0.3;

  // tau decay modes
  enum { hadron=1, electron, muon };

  // read input input file
  TChain chain("Delphes");
  chain.Add(inputfile);
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchGenJet = treeReader->UseBranch("GenJet");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  //TClonesArray *branchEvent = treeReader->UseBranch("Event");

  // set up loop variables
  GenParticle *genParticle;
  Jet *jet, *genJet;
  Electron *ele;
  Muon *mu;
  //LHEFEvent *event;

  // set up storage variables
  GenParticle *tau1, *tau2;
  Jet *genJetTau1, *genJetTau2;
  Jet *jetTau1, *jetTau2;
  Electron *eleTau1, *eleTau2;
  Muon *muTau1, *muTau2;
  GenParticle *b1, *b2;
  Jet *jetB1, *jetB2;

  // set up output variables and file
  Int_t nEvents;
  Float_t eventWeight;
  UInt_t tauDecayCat1, tauDecayCat2;
  UInt_t bTag1, bTag2;
  LorentzVector *genTau1=0, *genTau2=0, *genDecayTau1=0, *genDecayTau2=0, *recoTau1=0, *recoTau2=0;
  LorentzVector *genB1=0, *genB2=0, *recoB1=0, *recoB2=0;

  TFile *outFile = new TFile(outputfile, "RECREATE");

  TTree *sampTree = new TTree("Info", "Info");
  sampTree->Branch("nEvents",       &nEvents,        "nEvents/i");
  nEvents=numberOfEntries;
  sampTree->Fill();

  TTree *outTree = new TTree("Events", "Events");
  outTree->Branch("eventWeight",    &eventWeight,    "eventWeight/f");
  outTree->Branch("tauDecayCat1",   &tauDecayCat1,   "tauDecayCat1/i");
  outTree->Branch("tauDecayCat2",   &tauDecayCat2,   "tauDecayCat2/i");
  outTree->Branch("bTag1",          &bTag1,          "bTag1/i");
  outTree->Branch("bTag2",          &bTag2,          "bTag2/i");
  outTree->Branch("genTau1",        "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &genTau1);
  outTree->Branch("genTau2",        "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &genTau2);
  outTree->Branch("genDecayTau1",   "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &genDecayTau1);
  outTree->Branch("genDecayTau2",   "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &genDecayTau2);
  outTree->Branch("recoTau1",       "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &recoTau1);
  outTree->Branch("recoTau2",       "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &recoTau2);
  outTree->Branch("genB1",          "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &genB1);
  outTree->Branch("genB2",          "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &genB2);
  outTree->Branch("recoB1",         "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &recoB1);
  outTree->Branch("recoB2",         "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &recoB2);

  // set up more storage variables
  Int_t iGenTau1=-1,    iGenTau2=-1;
  Int_t iGenDecay1=-1,  iGenDecay2=-1;
  Int_t iRecoDecay1=-1, iRecoDecay2=-1;
  Int_t iGenB1=-1,      iGenB2=-1;
  Int_t iRecoB1=-1,     iRecoB2=-1;

  LorentzVector nothing(999,999,0,999);

  Int_t count=0;

  for (Int_t iEntry=0; iEntry<numberOfEntries; iEntry++) { // entry loop
    treeReader->ReadEntry(iEntry);

    //cout << " ---- " << endl;

    //event = (LHEFEvent*) branchEvent->At(0);
    //cout << xsec << ", " << numberOfEntries << ", " << event->Weight << ", ";
    eventWeight = 1;
    eventWeight = xsec;
    //eventWeight *= event->Weight;
    //cout << eventWeight << endl;

    iGenTau1=-1;    iGenTau2=-1;
    iGenDecay1=-1;  iGenDecay2=-1;
    iRecoDecay1=-1; iRecoDecay2=-1;
    iGenB1=-1;      iGenB2=-1;
    iRecoB1=-1;     iRecoB2=-1;

    tauDecayCat1=0; tauDecayCat2=0;
    bTag1=0;        bTag2=0;

    Int_t tCount=0, bCount=0;

    for (Int_t iParticle=0; iParticle<branchParticle->GetEntries(); iParticle++) {
      genParticle = (GenParticle*) branchParticle->At(iParticle);
      if ( fabs(genParticle->PID) == TAU_ID_CODE ) tCount++;
      if ( fabs(genParticle->PID) == B_ID_CODE) bCount++;
    }
    if ( ( tCount < 2 ) || ( bCount < 2 ) ) continue;
    //    cout << tCount << ", " << bCount << endl;

    // get generator level b and tau info
    for (Int_t iParticle=0; iParticle<branchParticle->GetEntries(); iParticle++) { // generator particle loop
      genParticle = (GenParticle*) branchParticle->At(iParticle);

      if ( fabs(genParticle->PID) == TAU_ID_CODE ) { // tau switch

	if (iGenTau1 == -1) {
	  iGenTau1 = iParticle;
	  tau1 = (GenParticle*) branchParticle->At(iGenTau1);
	}
	else if ( genParticle->PT > tau1->PT ) {
	  iGenTau2 = iGenTau1;
	  tau2 = tau1;
	  iGenTau1 = iParticle;
	  tau1 = (GenParticle*) branchParticle->At(iGenTau1);
	}
	else if ( iGenTau2 == -1 ) {
	  iGenTau2 = iParticle;
	  tau2 = (GenParticle*) branchParticle->At(iGenTau2);
	}
	else if ( genParticle->PT > tau2->PT) {
	  iGenTau2 = iParticle;
	  tau2 = (GenParticle*) branchParticle->At(iGenTau2);
	}
      } // end tau switch

      if ( fabs(genParticle->PID) == B_ID_CODE ) { // b switch

	if (iGenB1 == -1) {
	  iGenB1 = iParticle;
	  b1 = (GenParticle*) branchParticle->At(iGenB1);
	}
	else if ( genParticle->PT > b1->PT ) {
	  iGenB2 = iGenB1;
	  b2 = b1;
	  iGenB1 = iParticle;
	  b1 = (GenParticle*) branchParticle->At(iGenB2);
	}
	else if (iGenB2 == -1) {
	  iGenB2 = iParticle;
	  b2 = (GenParticle*) branchParticle->At(iGenB2);
	}
	else if ( genParticle->PT > b2->PT ) {
	  iGenB2 = iParticle;
	  b2 = (GenParticle*) branchParticle->At(iGenB2);	  
	}
      } // end b switch

    } // end generator particle loop

    // skip events without two generator level taus
    if ( (iGenTau1 == -1) || (iGenTau2 == -1) ) continue;

    // skip events without two generator level bs
    if ( (iGenB1 == -1) || (iGenB2 == -1) ) continue;
    
    // get generator level tau jet info
    for (Int_t iGenJet=0; iGenJet<branchGenJet->GetEntries(); iGenJet++) { // generator level jet loop
      genJet = (Jet*) branchGenJet->At(iGenJet);

      if ( deltaR(genJet->Eta, tau1->Eta, genJet->Phi, tau1->Phi) < MAX_MATCH_DIST ) {
	iGenDecay1 = iGenJet;
	genJetTau1 = (Jet*) branchGenJet->At(iGenDecay1);
      }
      else if ( deltaR(genJet->Eta, tau2->Eta, genJet->Phi, tau2->Phi) < MAX_MATCH_DIST ) {
	iGenDecay2 = iGenJet;
	genJetTau2 = (Jet*) branchGenJet->At(iGenDecay2);
      }
    } // end gen jet loop

    // get info for hadronic decay taus and b jets
    for (Int_t iJet=0; iJet<branchJet->GetEntries(); iJet++) { // reconstructed jet loop
      jet = (Jet*) branchJet->At(iJet);

      if ( jet->TauTag !=0 ) { // tau tag switch
	if ( deltaR(jet->Eta, tau1->Eta, jet->Phi, tau1->Phi) < MAX_MATCH_DIST ) {
	  //cout << "tau1 " << jet->TauTag << endl;
	  //cout << "(" << jet->Eta << " - " << tau1->Eta << ") ("  << jet->Phi << " - " << tau1->Phi << ")" << endl;
	  //cout << deltaR(jet->Eta, tau1->Eta, jet->Phi, tau1->Phi) << endl;
	  iRecoDecay1 = iJet;
	  tauDecayCat1 = hadron;
	  jetTau1 = (Jet*) branchJet->At(iRecoDecay1);
	}
	else if ( deltaR(jet->Eta, tau2->Eta, jet->Phi, tau2->Phi) < MAX_MATCH_DIST ) {
	  //cout << "tau2 " << jet->TauTag << endl;
	  //cout << "(" << jet->Eta << " - " << tau2->Eta << ") ("  << jet->Phi << " - " << tau2->Phi << ")" << endl;
	  //cout << deltaR(jet->Eta, tau2->Eta, jet->Phi, tau2->Phi) << endl;
	  iRecoDecay2 = iJet;
	  tauDecayCat2 = hadron;
	  jetTau2 = (Jet*) branchJet->At(iRecoDecay2);
	}
      } // end tau tag switch

      if ( jet->BTag !=0) { // b tag switch
	if ( deltaR(jet->Eta, b1->Eta, jet->Phi, b1->Phi) < MAX_MATCH_DIST ) {
	  //cout << "b1 " << jet->BTag << endl;
	  //cout << "(" << jet->Eta << " - " << b1->Eta << ") ("  << jet->Phi << " - " << b1->Phi << ")" << endl;
	  //cout << deltaR(jet->Eta, b1->Eta, jet->Phi, b1->Phi) << endl;
	  iRecoB1 = iJet;
	  bTag1 = jet->BTag;
	  jetB1 = (Jet*) branchJet->At(iRecoB1);
	}
	else if ( deltaR(jet->Eta, b2->Eta, jet->Phi, b2->Phi) < MAX_MATCH_DIST ) {
	  //cout << "b2 " << jet->BTag << endl;
	  //cout << "(" << jet->Eta << " - " << b2->Eta << ") ("  << jet->Phi << " - " << b2->Phi << ")" << endl;
	  //cout << deltaR(jet->Eta, b2->Eta, jet->Phi, b2->Phi) << endl;
	  iRecoB2 = iJet;
	  bTag2 = jet->BTag;
	  jetB2 = (Jet*) branchJet->At(iRecoB2);
	}
      } // end b tag switch

    } // end reco jet loop

    // get info for taus that decay to muons
    for (Int_t iMuon=0; iMuon<branchMuon->GetEntries(); iMuon++) { // reco muon loop
      mu = (Muon*) branchMuon->At(iMuon);

      if ( deltaR(mu->Eta, tau1->Eta, mu->Phi, tau1->Phi) < MAX_MATCH_DIST ) {
	if ( iRecoDecay1!=-1 ) { cout << "found a muon and a jet?" << endl; continue; }
	iRecoDecay1 = iMuon;
	tauDecayCat1 = muon;
	muTau1 = (Muon*) branchMuon->At(iRecoDecay1);
      }
      else if ( deltaR(mu->Eta, tau2->Eta, mu->Phi, tau2->Phi) < MAX_MATCH_DIST ) {
	if ( iRecoDecay2!=-1 ) { cout << "found a muon and a jet?" << endl; continue; }
	iRecoDecay2 = iMuon;
	tauDecayCat2 = muon;
	muTau2 = (Muon*) branchMuon->At(iRecoDecay2);
      }
    } // end muon loop

    // get info for taus that decay to electrons
    for (Int_t iEle=0; iEle<branchElectron->GetEntries(); iEle++) { // reco ele loop
      ele = (Electron*) branchElectron->At(iEle);

      if ( deltaR(ele->Eta, tau1->Eta, ele->Phi, tau1->Phi) < MAX_MATCH_DIST ) {
	if ( iRecoDecay1!=-1) { cout << "found an electron and something else?" << endl; continue; }
	iRecoDecay1 = iEle;
	tauDecayCat1 = electron;
	eleTau1 = (Electron*) branchElectron->At(iRecoDecay1);
      }
      else if ( deltaR(ele->Eta, tau2->Eta, ele->Phi, tau2->Phi) < MAX_MATCH_DIST ) {
	if ( iRecoDecay2!=-1) { cout << "found an electron and something else?" << endl; continue; }
	iRecoDecay2 = iEle;
	tauDecayCat2 = electron;
	eleTau2 = (Electron*) branchElectron->At(iRecoDecay2);
      }
    } // end electron loop

    // 
    // OUTPUT 
    //

    // store generator particles
    LorentzVector vGenTau1(tau1->PT, tau1->Eta, tau1->Phi, tau1->Mass);
    genTau1 = &vGenTau1;
    LorentzVector vGenTau2(tau2->PT, tau2->Eta, tau2->Phi, tau2->Mass);
    genTau2 = &vGenTau2;

    LorentzVector vGenB1(b1->PT, b1->Eta, b1->Phi, b1->Mass);
    genB1 = &vGenB1;
    LorentzVector vGenB2(b2->PT, b2->Eta, b2->Phi, b2->Mass);
    genB2 = &vGenB2;

    // store generator jets
    LorentzVector vGenJetTau1(0,0,0,0);
    if ( iGenDecay1 == -1 ) genDecayTau1 = &nothing;
    else {
      vGenJetTau1.SetPt(genJetTau1->PT);
      vGenJetTau1.SetEta(genJetTau1->Eta);
      vGenJetTau1.SetPhi(genJetTau1->Phi);
      vGenJetTau1.SetM(genJetTau1->Mass);
      genDecayTau1 = &vGenJetTau1;
    }
    LorentzVector vGenJetTau2(0,0,0,0);
    if ( iGenDecay2 == -1 ) genDecayTau2 = &nothing;
    else {
      vGenJetTau2.SetPt(genJetTau2->PT);
      vGenJetTau2.SetEta(genJetTau2->Eta);
      vGenJetTau2.SetPhi(genJetTau2->Phi);
      vGenJetTau2.SetM(genJetTau2->Mass);
      genDecayTau2 = &vGenJetTau2;
    }

    LorentzVector vRecoB1(0,0,0,0);
    if (bTag1==0) recoB1 = &nothing;
    else {
      vRecoB1.SetPt(jetB1->PT);
      vRecoB1.SetEta(jetB1->Eta);
      vRecoB1.SetPhi(jetB1->Phi);
      vRecoB1.SetM(jetB1->Mass);
      recoB1 = &vRecoB1;
    }

    LorentzVector vRecoB2(0,0,0,0);
    if (bTag2==0) recoB2 = &nothing;
    else {
      vRecoB2.SetPt(jetB2->PT);
      vRecoB2.SetEta(jetB2->Eta);
      vRecoB2.SetPhi(jetB2->Phi);
      vRecoB2.SetM(jetB2->Mass);
      recoB2 = &vRecoB2;
    }

    LorentzVector vRecoTau1(0,0,0,0);
    if ( tauDecayCat1 == 0 ) recoTau1 = &nothing;
    else if ( tauDecayCat1 == hadron ) {
      vRecoTau1.SetPt(jetTau1->PT);
      vRecoTau1.SetEta(jetTau1->Eta);
      vRecoTau1.SetPhi(jetTau1->Phi);
      vRecoTau1.SetM(jetTau1->Mass);
      recoTau1 = &vRecoTau1;
    }
    else if ( tauDecayCat1 == muon ) {
      vRecoTau1.SetPt(muTau1->PT);
      vRecoTau1.SetEta(muTau1->Eta);
      vRecoTau1.SetPhi(muTau1->Phi);
      vRecoTau1.SetM(MUON_MASS);
      recoTau1 = &vRecoTau1;
    }
    else if ( tauDecayCat1 == electron ) {
      vRecoTau1.SetPt(eleTau1->PT);
      vRecoTau1.SetEta(eleTau1->Eta);
      vRecoTau1.SetPhi(eleTau1->Phi);
      vRecoTau1.SetM(ELE_MASS);
      recoTau1 = &vRecoTau1;
    }

    LorentzVector vRecoTau2(0,0,0,0);
    if ( tauDecayCat2 == 0 ) recoTau2 = &nothing;
    else if ( tauDecayCat2 == hadron ) {
      vRecoTau2.SetPt(jetTau2->PT);
      vRecoTau2.SetEta(jetTau2->Eta);
      vRecoTau2.SetPhi(jetTau2->Phi);
      vRecoTau2.SetM(jetTau2->Mass);
      recoTau2 = &vRecoTau2;
    }
    else if ( tauDecayCat2 == muon ) {
      vRecoTau2.SetPt(muTau2->PT);
      vRecoTau2.SetEta(muTau2->Eta);
      vRecoTau2.SetPhi(muTau2->Phi);
      vRecoTau2.SetM(MUON_MASS);
      recoTau2 = &vRecoTau2;
    }
    else if ( tauDecayCat2 == electron ) {
      vRecoTau2.SetPt(eleTau2->PT);
      vRecoTau2.SetEta(eleTau2->Eta);
      vRecoTau2.SetPhi(eleTau2->Phi);
      vRecoTau2.SetM(ELE_MASS);
      recoTau2 = &vRecoTau2;
    }

    outTree->Fill();

    if ( ( bTag1==0 ) || (bTag2==0) || (tauDecayCat1==0) || (tauDecayCat2==0) ) continue;

    count++;
    cout << count << "/" << nEvents << endl;

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
