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

void select_vbf(const TString inputfile="root://eoscms.cern.ch//store/group/phys_higgs/upgrade/PhaseII/Configuration4v2/140PileUp/Bjj-vbf-4p-0-700-v1510_14TEV/Bjj-vbf-4p-0-700-v1510_14TEV_265719429_PhaseII_Conf4v2_140PileUp.root",
	       const Float_t xsec=1.0,
	       const TString outputfile="test.root") {

  // declare constants
  const Double_t MUON_MASS = 0.105658369;
  const Double_t ELE_MASS  = 0.000511;

  const Int_t TAU_ID_CODE = 15;
  const Int_t B_ID_CODE = 5;

  const Int_t T_ID_CODE = 6;
  const Int_t Z_ID_CODE = 23;
  const Int_t W_ID_CODE = 24;
  const Int_t H_ID_CODE = 25;

  const Float_t MAX_MATCH_DIST = 0.5;

  // event types
  enum { HH=0, TT, ZH, WH, WW, ZZ, ZW, ETC };

  // tau decay modes
  enum { hadron=1, electron, muon };

  // read input input file
  TChain chain("Delphes");
  chain.Add(inputfile);
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchMET =treeReader->UseBranch("MissingET");

  TClonesArray *branchGenJet = treeReader->UseBranch("GenJet");
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchEvent = treeReader->UseBranch("Event");

  //set up loop variables
  GenParticle *genParticle=0;
  Jet *genJet=0;
  Jet *jet=0;
  Electron *ele=0;
  Muon *mu=0;
  MissingET *missET=0;
  LHEFEvent *event=0;

  // set up storage variables
  Jet *jetTau1=0, *jetTau2=0;
  Electron *eleTau=0;
  Muon *muTau=0;
  Jet *jetB1=0, *jetB2=0;
  GenParticle *genTau1=0, *genTau2=0;
  GenParticle *genB1=0, *genB2=0;
  Jet *genJetB1=0, *genJetB2=0;
  Jet *genJetTau1=0, *genJetTau2=0;

  Int_t iB1=-1,         iB2=-1;
  Int_t iT1=-1,         iT2=-1;
  Int_t iGenTau1=-1,    iGenTau2=-1;
  Int_t iGenB1=-1,      iGenB2=-1;
  Int_t iGenJetTau1=-1, iGenJetTau2=-1;
  Int_t iGenJetB1=-1,   iGenJetB2=-1;

  // set up output variables and file
  Int_t nEvents;
  Float_t eventWeight;

  Float_t met, metPhi;

  Int_t eventType;
  UInt_t tauCat1=0, tauCat2=0;
  UInt_t bTag1=0, bTag2=0;

  LorentzVector *sRecoTau1=0, *sRecoTau2=0;
  LorentzVector *sGenJetTau1=0, *sGenJetTau2=0;
  LorentzVector *sGenTau1=0, *sGenTau2=0;

  LorentzVector *sRecoB1=0, *sRecoB2=0;
  LorentzVector *sGenJetB1=0, *sGenJetB2=0;
  LorentzVector *sGenB1=0, *sGenB2=0;

  TFile *outFile = new TFile(outputfile, "RECREATE");

  // tree to hold the number of events in the file before selection
  TTree *sampTree = new TTree("Info", "Info");
  sampTree->Branch("nEvents",       &nEvents,        "nEvents/i");
  nEvents=numberOfEntries;
  sampTree->Fill();

  // tree to hold information about selected events
  TTree *outTree = new TTree("Events", "Events");
  outTree->Branch("eventWeight",    &eventWeight,    "eventWeight/f");  // event weight from cross-section and Event->Weight
  outTree->Branch("eventType",      &eventType,      "eventType/i");    // event type (0=signal, 1=tt, 2=zh, 3=other)
  outTree->Branch("tauCat1",        &tauCat1,        "tauCat1/i");      // leading tau final state - jet, muon, electron
  outTree->Branch("tauCat2",        &tauCat2,        "tauCat2/i");      // second tau final state - jet, muon, electron
  outTree->Branch("bTag1",          &bTag1,          "bTag1/i");        // leading b-jet tag from delphes
  outTree->Branch("bTag2",          &bTag2,          "bTag2/i");        // second b-jet tag from delphes
  outTree->Branch("met",            &met,            "met/f");          // missing transverse energy
  outTree->Branch("metPhi",         &metPhi,         "metPhi/f");       // missing transverse energy phi
  outTree->Branch("sGenTau1",       "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sGenTau1);      // 4-vector for generator leading tau
  outTree->Branch("sGenTau2",       "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sGenTau2);      // 4-vector for generator second tau
  outTree->Branch("sGenB1",         "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sGenB1);        // 4-vector for generator leading b-jet
  outTree->Branch("sGenB2",         "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sGenB2);        // 4-vector for generator second b-jet
  outTree->Branch("sGenJetTau1",    "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sGenJetTau1);   // 4-vector for generator leading tau
  outTree->Branch("sGenJetTau2",    "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sGenJetTau2);   // 4-vector for generator second tau
  outTree->Branch("sGenJetB1",      "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sGenJetB1);     // 4-vector for generator leading b-jet
  outTree->Branch("sGenJetB2",      "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sGenJetB2);     // 4-vector for generator second b-jet
  outTree->Branch("sRecoTau1",      "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sRecoTau1);     // 4-vector for reconstructed leading tau
  outTree->Branch("sRecoTau2",      "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sRecoTau2);     // 4-vector for reconstructed second tau
  outTree->Branch("sRecoB1",        "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sRecoB1);       // 4-vector for reconstructed leading b-jet
  outTree->Branch("sRecoB2",        "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sRecoB2);       // 4-vector for reconstructed second b-jet

  // define placeholder vector for things that don't exist
  LorentzVector nothing(999,999,0,999);

  Float_t select=0;

  for (Int_t iEntry=0; iEntry<numberOfEntries; iEntry++) { // entry loop
    treeReader->ReadEntry(iEntry);

    // ********************
    // RESET
    // ********************
    
    iB1=-1; iB2=-1; iT1=-1; iT2=-1;
    iGenTau1=-1;    iGenTau2=-1;
    iGenB1=-1;      iGenB2=-1;
    iGenJetTau1=-1; iGenJetTau2=-1;
    iGenJetB1=-1;   iGenJetB2=-1;
    tauCat1=-1; tauCat2=-1; bTag1=-1; bTag2=-1;
    eventType=-1;

    jetTau1=0; jetTau2=0; eleTau=0; muTau=0;
    jetB1=0;   jetB2=0;
    genTau1=0; genTau2=0; genB1=0;  genB2=0;
    genJetB1=0; genJetB2=0; genJetTau1=0; genJetTau2=0;
    sGenTau1=0; sGenTau2=0; sGenB1=0;  sGenB2=0;
    sGenJetB1=0; sGenJetB2=0; sGenJetTau1=0; sGenJetTau2=0;

    // ********************
    // EVENT WEIGHT
    // ********************

    eventWeight = 1;
    if (branchEvent) {
      event = (LHEFEvent*) branchEvent->At(0);
      eventWeight*=event->Weight;
    }
    eventWeight *= xsec;

    // ********************
    // RECO OBJECTS
    // ********************

    // get reconstructed hadronic taus
    for (Int_t iJet=0; iJet<branchJet->GetEntries(); iJet++) { // reconstructed jet loop
      jet = (Jet*) branchJet->At(iJet);

      if (jet->PT < 30) continue;
      if (fabs(jet->Eta)>2.5) continue;
      if (jet->TauTag==0) continue;
      if ((jetTau1)&&(deltaR(jet->Eta, jetTau1->Eta, jet->Phi, jetTau1->Phi) < MAX_MATCH_DIST)) continue;
      if ((jetTau2)&&(deltaR(jet->Eta, jetTau2->Eta, jet->Phi, jetTau2->Phi) < MAX_MATCH_DIST)) continue;

      if (iT1==-1) { 
	iT1=iJet; 
	jetTau1 = (Jet*) branchJet->At(iT1); 
	tauCat1=hadron;
      }
      else if (jet->PT > jetTau1->PT) { 
	iT2=iT1; 
	jetTau2 = (Jet*) branchJet->At(iT2); 
	tauCat2=hadron; 
	iT1=iJet; 
	jetTau1 = (Jet*) branchJet->At(iT1); 
	tauCat1=hadron;
      }
      else if (iT2==-1) { 
	iT2=iJet; 
	jetTau2 = (Jet*) branchJet->At(iT2); 
	tauCat2=hadron;
      }
      else if (jet->PT > jetTau2->PT) { 
	iT2=iJet; 
	jetTau2 = (Jet*) branchJet->At(iT2); 
	tauCat2=hadron;
      }
    } // end reco jet loop

    for (Int_t iJet=0; iJet<branchJet->GetEntries(); iJet++) { // reconstructed jet loop
      jet = (Jet*) branchJet->At(iJet);

      if (jet->PT < 30) continue;
      if (fabs(jet->Eta)>4.7) continue;
      //if (jet->BTag==0) continue;
      if ((jetTau1)&&(deltaR(jet->Eta, jetTau1->Eta, jet->Phi, jetTau1->Phi) < MAX_MATCH_DIST)) continue;
      if ((jetTau2)&&(deltaR(jet->Eta, jetTau2->Eta, jet->Phi, jetTau2->Phi) < MAX_MATCH_DIST)) continue;

      if (iB1==-1) {
	iB1=iJet; 
	jetB1 = (Jet*) branchJet->At(iB1); 
	bTag1=jetB1->BTag;
      }
      else if (jet->PT > jetB1->PT) {
	iB2=iB1; 
	jetB2 = (Jet*) branchJet->At(iB2); 
	bTag2=bTag1; 
	iB1=iJet;
	jetB1 = (Jet*) branchJet->At(iB1);
	bTag1=jetB1->BTag;
      }
      else if (iB2==-1) { 
	iB2=iJet; 
	jetB2 = (Jet*) branchJet->At(iB2); 
	bTag2=jetB2->BTag;
      }
      else if (jet->PT > jetB2->PT) { 
	iB2=iJet; 
	jetB2 = (Jet*) branchJet->At(iB1); 
	bTag2=jetB2->BTag;
      }
    }
    
    if ((iT1==-1) || (iT2==-1)) {
      // get muonic taus
      for (Int_t iMuon=0; iMuon<branchMuon->GetEntries(); iMuon++) { // reco muon loop
	mu = (Muon*) branchMuon->At(iMuon);
	if (fabs(mu->Eta)>2.5) continue;
	if (mu->PT < 30) continue;

	if ((jetTau1)&&(deltaR(mu->Eta, jetTau1->Eta, mu->Phi, jetTau1->Phi) < MAX_MATCH_DIST)) continue;
	if ((jetTau2)&&(deltaR(mu->Eta, jetTau2->Eta, mu->Phi, jetTau2->Phi) < MAX_MATCH_DIST)) continue;

	if (iT1==-1) { 
	  iT1=iMuon; 
	  muTau = (Muon*) branchMuon->At(iT1); 
	  tauCat1=muon; 
	}
	else if (iT2==-1) { 
	  iT2=iMuon; 
	  muTau = (Muon*) branchMuon->At(iT2); 
	  tauCat2=muon; 
	}
	else if (muTau) {
	  if ( mu->PT > muTau->PT ) { 
	    if (tauCat1==muon) {
	      iT1=iMuon;
	      muTau = (Muon*) branchMuon->At(iT1); 
	    }
	    else if (tauCat2==muon) {
	      iT2=iMuon;
	      muTau = (Muon*) branchMuon->At(iT2); 
	    }
	  }
	}
      }

      // get electronic taus
      for (Int_t iEle=0; iEle<branchElectron->GetEntries(); iEle++) { // reco ele loop
	ele = (Electron*) branchElectron->At(iEle);

	if (ele->PT < 30) continue;
	if (fabs(ele->Eta)>2.5) continue;
	if ((jetTau1)&&(deltaR(ele->Eta, jetTau1->Eta, ele->Phi, jetTau1->Phi) < MAX_MATCH_DIST)) continue;
	if ((jetTau2)&&(deltaR(ele->Eta, jetTau2->Eta, ele->Phi, jetTau2->Phi) < MAX_MATCH_DIST)) continue;
	
	if (iT1==-1) { 
	  iT1=iEle; 
	  eleTau = (Electron*) branchElectron->At(iT1); 
	  tauCat1=electron; 
	}
	else if (iT2==-1) { 
	  iT2=iEle; 
	  eleTau = (Electron*) branchElectron->At(iT2); 
	  tauCat2=electron;
	}
	else if (eleTau) {
	  if ( ele->PT > eleTau->PT ) { 
	    if ( tauCat1==electron) {
	      iT1=iEle;
	      eleTau = (Electron*) branchElectron->At(iT1); 
	    }
	    else if ( tauCat2==electron) {
	      iT2=iEle;
	      eleTau = (Electron*) branchElectron->At(iT2); 
	    }
	  }
	}
      }
    }
      
    if ((iT1==-1)||(iT2==-1)||(iB1==-1)||(iB2==-1)) continue;

    // fill 4-vector for leading b-jet
    LorentzVector vRecoB1(0,0,0,0);
    if (jetB1) {
      vRecoB1.SetPt(jetB1->PT);
      vRecoB1.SetEta(jetB1->Eta);
      vRecoB1.SetPhi(jetB1->Phi);
      vRecoB1.SetM(jetB1->Mass);
      sRecoB1 = &vRecoB1;
    }
    else sRecoB1 = &nothing;
    // fill 4-vector for second b-jet
    LorentzVector vRecoB2(0,0,0,0);
    if (jetB2) {
      vRecoB2.SetPt(jetB2->PT);
      vRecoB2.SetEta(jetB2->Eta);
      vRecoB2.SetPhi(jetB2->Phi);
      vRecoB2.SetM(jetB2->Mass);
      sRecoB2 = &vRecoB2;
    }
    else sRecoB2 = &nothing;

    LorentzVector dijet=vRecoB1+vRecoB2;
    if (dijet.M() < 500) continue;

    if (fabs(vRecoB1.Eta()-vRecoB2.Eta())<3.5) continue;

    Int_t centJet=0;

    for (Int_t iJet=0; iJet<branchJet->GetEntries(); iJet++) { // reconstructed jet loop
      jet = (Jet*) branchJet->At(iJet);

      if (jet->PT < 30) continue;
      if (fabs(jet->Eta)>4.7) continue;
      if ((iJet==iB1)||(iJet==iB2)) continue;
      if ((jetTau1)&&(deltaR(mu->Eta, jetTau1->Eta, mu->Phi, jetTau1->Phi) < MAX_MATCH_DIST)) continue;
      if ((jetTau2)&&(deltaR(mu->Eta, jetTau2->Eta, mu->Phi, jetTau2->Phi) < MAX_MATCH_DIST)) continue;

      if ((jetB1->Eta > jetB2->Eta) && (jet->Eta > jetB2->Eta) && (jet->Eta < jetB2->Eta)) centJet++;
      else if ((jetB2->Eta > jetB1->Eta) && (jet->Eta > jetB1->Eta) && (jet->Eta < jetB2->Eta)) centJet++;
    }

    if (centJet!=0) continue;

    // fill 4-vector for leading tau
    LorentzVector vRecoTau1(0,0,0,0);
    if (jetTau1) {
      vRecoTau1.SetPt(jetTau1->PT);
      vRecoTau1.SetEta(jetTau1->Eta);
      vRecoTau1.SetPhi(jetTau1->Phi);
      vRecoTau1.SetM(jetTau1->Mass);
      sRecoTau1 = &vRecoTau1;
    }
    else if ((muTau)&&(tauCat1==muon)) {
      vRecoTau1.SetPt(muTau->PT);
      vRecoTau1.SetEta(muTau->Eta);
      vRecoTau1.SetPhi(muTau->Phi);
      vRecoTau1.SetM(MUON_MASS);
      sRecoTau1 = &vRecoTau1;
    }
    else if ((eleTau)&&(tauCat1==electron)) {
      vRecoTau1.SetPt(eleTau->PT);
      vRecoTau1.SetEta(eleTau->Eta);
      vRecoTau1.SetPhi(eleTau->Phi);
      vRecoTau1.SetM(ELE_MASS);
      sRecoTau1 = &vRecoTau1;
    }
    else sRecoTau1 = &nothing;

    // fill 4-vector for second tau
    LorentzVector vRecoTau2(0,0,0,0);
    if (jetTau2) {
      vRecoTau2.SetPt(jetTau2->PT);
      vRecoTau2.SetEta(jetTau2->Eta);
      vRecoTau2.SetPhi(jetTau2->Phi);
      vRecoTau2.SetM(jetTau2->Mass);
      sRecoTau2 = &vRecoTau2;
    }
    else if ((muTau)&&(tauCat2==muon)) {
      vRecoTau2.SetPt(muTau->PT);
      vRecoTau2.SetEta(muTau->Eta);
      vRecoTau2.SetPhi(muTau->Phi);
      vRecoTau2.SetM(MUON_MASS);
      sRecoTau2 = &vRecoTau2;
    }
    else if ((eleTau)&&(tauCat2==electron)) {
      vRecoTau2.SetPt(eleTau->PT);
      vRecoTau2.SetEta(eleTau->Eta);
      vRecoTau2.SetPhi(eleTau->Phi);
      vRecoTau2.SetM(ELE_MASS);
      sRecoTau2 = &vRecoTau2;
    }
    else sRecoTau2 = &nothing;

    // ********************
    // GEN PARTICLES
    // ********************

    for (Int_t iParticle=0; iParticle<branchParticle->GetEntries(); iParticle++) { // generator particle loop
      genParticle = (GenParticle*) branchParticle->At(iParticle);

      if ( fabs(genParticle->PID) == TAU_ID_CODE ) { // tau switch
	if ( (tauCat1 != -1) && ( deltaR(genParticle->Eta, vRecoTau1.Eta(), genParticle->Phi, vRecoTau1.Phi()) < MAX_MATCH_DIST )) {
	  iGenTau1 = iParticle;
	  genTau1 = (GenParticle*) branchParticle->At(iGenTau1);
	}
	else if ( (tauCat2 != -1) && ( deltaR(genParticle->Eta, vRecoTau2.Eta(), genParticle->Phi, vRecoTau2.Phi()) < MAX_MATCH_DIST )) {
	  iGenTau2 = iParticle;
	  genTau2 = (GenParticle*) branchParticle->At(iGenTau2);
	}
	else if ( (tauCat1 == -1) && (iGenTau1==-1) ) { 
	  iGenTau1 = iParticle;
	  genTau1 = (GenParticle*) branchParticle->At(iGenTau1);
	}
	else if ( (tauCat2 == -1) && (iGenTau1==-2) ) { 
	  iGenTau2 = iParticle;
	  genTau2 = (GenParticle*) branchParticle->At(iGenTau2);
	}
      }
      /*
      if ( fabs(genParticle->PID) == B_ID_CODE ) { // b switch

	if ( (bTag1 != -1) && ( deltaR(genParticle->Eta, vRecoB1.Eta(), genParticle->Phi, vRecoB1.Phi()) < MAX_MATCH_DIST )) {
	  iGenB1 = iParticle;
	  genB1 = (GenParticle*) branchParticle->At(iGenB1);
	}
	else if ( (bTag2 != -1) && ( deltaR(genParticle->Eta, vRecoB2.Eta(), genParticle->Phi, vRecoB2.Phi()) < MAX_MATCH_DIST )) {
	  iGenB2 = iParticle;
	  genB2 = (GenParticle*) branchParticle->At(iGenB2);
	}
	else if ( (bTag1 == -1) && (iGenB1==-1) ) { 
	  iGenB1 = iParticle;
	  genB1 = (GenParticle*) branchParticle->At(iGenB1);
	}
	else if ( (bTag2 == -1) && (iGenB1==-2) ) { 
	  iGenB2 = iParticle;
	  genB2 = (GenParticle*) branchParticle->At(iGenB2);
	}
      }
      */
    }

    LorentzVector vGenTau1(0,0,0,0);
    if (genTau1) {
      vGenTau1.SetPt(genTau1->PT);
      vGenTau1.SetEta(genTau1->Eta);
      vGenTau1.SetPhi(genTau1->Phi);
      vGenTau1.SetM(genTau1->Mass);
      sGenTau1 = &vGenTau1;
    }
    else sGenTau1 = &nothing;

    LorentzVector vGenTau2(0,0,0,0);
    if (genTau2) {
      vGenTau2.SetPt(genTau2->PT);
      vGenTau2.SetEta(genTau2->Eta);
      vGenTau2.SetPhi(genTau2->Phi);
      vGenTau2.SetM(genTau2->Mass);
      sGenTau2 = &vGenTau2;
    }
    else sGenTau2 = &nothing;

    LorentzVector vGenB1(0,0,0,0);
    if (genB1) {
      vGenB1.SetPt(genB1->PT);
      vGenB1.SetEta(genB1->Eta);
      vGenB1.SetPhi(genB1->Phi);
      vGenB1.SetM(genB1->Mass);
      sGenB1 = &vGenB1;
    }
    else sGenB1 = &nothing;

    LorentzVector vGenB2(0,0,0,0);
    if (genB2) {
      vGenB2.SetPt(genB2->PT);
      vGenB2.SetEta(genB2->Eta);
      vGenB2.SetPhi(genB2->Phi);
      vGenB2.SetM(genB2->Mass);
      sGenB2 = &vGenB2;
    }
    else sGenB2 = &nothing;

    // match generator level jets to generator particles
    for (Int_t iJet=0; iJet<branchGenJet->GetEntries(); iJet++) { // generator level jet loop
      genJet = (Jet*) branchGenJet->At(iJet);

      if ((genB1) && (deltaR(genJet->Eta, genB1->Eta, genJet->Phi, genB1->Phi) < MAX_MATCH_DIST) ) {
	iGenJetB1=iJet;
	genJetB1 = (Jet*) branchGenJet->At(iGenJetB1);
      }

      else if ((genB2) && (deltaR(genJet->Eta, genB2->Eta, genJet->Phi, genB2->Phi) < MAX_MATCH_DIST) ) {
	iGenJetB2=iJet;
	genJetB2 = (Jet*) branchGenJet->At(iGenJetB2);
      }

      else if ((genTau1) && (deltaR(genJet->Eta, genTau1->Eta, genJet->Phi, genTau1->Phi) < MAX_MATCH_DIST) ) {
	iGenJetTau1=iJet;
	genJetTau1 = (Jet*) branchGenJet->At(iGenJetTau1);
      }

      else if ((genTau2) && (deltaR(genJet->Eta, genTau2->Eta, genJet->Phi, genTau2->Phi) < MAX_MATCH_DIST) ) {
	iGenJetTau2=iJet;
	genJetTau2 = (Jet*) branchGenJet->At(iGenJetTau2);
      }

    }

    LorentzVector vGenJetTau1(0,0,0,0);
    if (genJetTau1) {
      vGenJetTau1.SetPt(genJetTau1->PT);
      vGenJetTau1.SetEta(genJetTau1->Eta);
      vGenJetTau1.SetPhi(genJetTau1->Phi);
      vGenJetTau1.SetM(genJetTau1->Mass);
      sGenJetTau1 = &vGenJetTau1;
    }
    else sGenJetTau1 = &nothing;

    LorentzVector vGenJetTau2(0,0,0,0);
    if (genJetTau2) {
      vGenJetTau2.SetPt(genJetTau2->PT);
      vGenJetTau2.SetEta(genJetTau2->Eta);
      vGenJetTau2.SetPhi(genJetTau2->Phi);
      vGenJetTau2.SetM(genJetTau2->Mass);
      sGenJetTau2 = &vGenJetTau2;
    }
    else sGenJetTau2 = &nothing;

    LorentzVector vGenJetB1(0,0,0,0);
    if (genJetB1) {
      vGenJetB1.SetPt(genJetB1->PT);
      vGenJetB1.SetEta(genJetB1->Eta);
      vGenJetB1.SetPhi(genJetB1->Phi);
      vGenJetB1.SetM(genJetB1->Mass);
      sGenJetB1 = &vGenJetB1;
    }
    else sGenJetB1 = &nothing;

    LorentzVector vGenJetB2(0,0,0,0);
    if (genJetB2) {
      vGenJetB2.SetPt(genJetB2->PT);
      vGenJetB2.SetEta(genJetB2->Eta);
      vGenJetB2.SetPhi(genJetB2->Phi);
      vGenJetB2.SetM(genJetB2->Mass);
      sGenJetB2 = &vGenJetB2;
    }
    else sGenJetB2 = &nothing;

    // ********************
    // MT2 CALC
    // ********************

    missET = (MissingET*) branchMET->At(0);

    met=missET->MET;
    metPhi=missET->Phi;

    // ********************
    // BKGD SORTING
    // ********************

    Int_t nH=0, nW=0, nZ=0, nT=0;
    for (Int_t iPart=0; iPart<branchParticle->GetEntries(); iPart++) {
      genParticle = (GenParticle*) branchParticle->At(iPart);
      if (fabs(genParticle->PID)==H_ID_CODE) nH++;
      else if (fabs(genParticle->PID)==Z_ID_CODE) nZ++;
      else if (fabs(genParticle->PID)==W_ID_CODE) nW++;
      else if (fabs(genParticle->PID)==T_ID_CODE) nT++;
    }

    if ( (nH==2) && (nW==0) && (nZ==0) ) eventType=HH;
    else if ( (nH==0) && (nT==2) ) eventType=TT;
    else if ( (nH!=0) && (nZ>0) ) eventType=ZH;
    else if ( (nH!=0) && (nW>0) && (nZ==0) ) eventType=WH;
    else if ( (nH==0) && (nW>0) && (nZ>0) ) eventType=ZW;
    else if (nZ>1) eventType=ZZ;
    else if (nW>1) eventType=WW;
    else eventType=ETC;

    select+=eventWeight;

    outTree->Fill();

  } // end event loop

  cout << select << " of " << numberOfEntries << endl;

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
