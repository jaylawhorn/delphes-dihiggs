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

void selection(const TString inputfile="root://eoscms.cern.ch//store/group/phys_higgs/upgrade/PhaseII/Configuration4v2/140PileUp/HHToTTBB_14TeV/HHToTTBB_14TeV_1.root", 
	       const Float_t xsec=2.92,
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

  const Float_t MAX_MATCH_DIST = 0.3;

  // event types
  enum { HH=0, TT, ZH, WH, WW, ZZ, ZW, ETC };

  // tau decay modes
  enum { hadron=1, electron, muon };

  // setup mt2 minimizer
  ROOT::Math::Minimizer *min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
  min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
  min->SetTolerance(10.0);
  min->SetPrintLevel(0);

  TVector2 tau1(0,0), tau2(0,0), mpt(0,0);
  TVector2 b1(0,0), b2(0,0);
  Float_t mTau1=0, mTau2=0;
  Float_t mB1=0, mB2=0;
  Double_t mt2=0;

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
  //TClonesArray *branchEvent = treeReader->UseBranch("Event");

  //set up loop variables
  GenParticle *genParticle;
  Jet *genJet;
  Jet *jet;
  Electron *ele;
  Muon *mu;
  MissingET *missET;
  //LHEFEvent *event;

  // set up storage variables
  Jet *jetTau1, *jetTau2;
  Electron *eleTau;
  Muon *muTau;
  Jet *jetB1, *jetB2;
  Jet *extraJet;
  GenParticle *genTau1, *genTau2;
  GenParticle *genB1, *genB2;
  Jet *genJetB1, *genJetB2;
  Jet *genJetTau1, *genJetTau2;

  Int_t iB1=-1,         iB2=-1;
  Int_t iT1=-1,         iT2=-1;
  Int_t iGenTau1=-1,    iGenTau2=-1;
  Int_t iGenB1=-1,      iGenB2=-1;
  Int_t iGenJetTau1=-1, iGenJetTau2=-1;
  Int_t iGenJetB1=-1,   iGenJetB2=-1;
  Int_t iE=-1;

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

  LorentzVector *sRecoExtra=0;

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
  outTree->Branch("mt2",            &mt2,            "mt2/D");          // "stransverse mass"
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
  outTree->Branch("sRecoExtra",     "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sRecoExtra);    // 4-vector for reconstructed extra jet

  // define placeholder vector for things that don't exist
  LorentzVector nothing(999,999,0,999);

  for (Int_t iEntry=0; iEntry<numberOfEntries; iEntry++) { // entry loop
    treeReader->ReadEntry(iEntry);

    // ********************
    // RESET
    // ********************

    iB1=-1; iB2=-1; iT1=-1; iT2=-1; iE=-1;
    iGenTau1=-1;    iGenTau2=-1;
    iGenB1=-1;      iGenB2=-1;
    iGenJetTau1=-1; iGenJetTau2=-1;
    iGenJetB1=-1;   iGenJetB2=-1;

    tauCat1=0; tauCat2=0; bTag1=0; bTag2=0;
    eventType=-1;

    // ********************
    // EVENT WEIGHT
    // ********************

    // comment out following line for di-higgs samples
    //event = (LHEFEvent*) branchEvent->At(0);
    eventWeight = 1;
    eventWeight *= xsec;
    // comment out following line for di-higgs samples
    //eventWeight *= event->Weight;

    // ********************
    // RECO OBJECTS
    // ********************

    // get reconstructed hadronic taus and b jets
    for (Int_t iJet=0; iJet<branchJet->GetEntries(); iJet++) { // reconstructed jet loop
      jet = (Jet*) branchJet->At(iJet);

      if (jet->BTag!=0) {

	if ((iB1==-1)) { jetB1 = (Jet*) branchJet->At(iJet); iB1=iJet; bTag1=jetB1->BTag; }
	else if ((iB1!=iJet) && (jet->PT > jetB1->PT)) { jetB2 = (Jet*) branchJet->At(iB1); iB2=iB1; bTag2=bTag1; jetB1 = (Jet*) branchJet->At(iJet); iB1=iJet; bTag1=jetB1->BTag; }
	else if ((iB1!=iJet) && (iB2==-1)) { jetB2 = (Jet*) branchJet->At(iJet); iB2=iJet; bTag2=jetB2->BTag; }
	else if ((iB1!=iJet) && (jet->PT > jetB2->PT)) { jetB2 = (Jet*) branchJet->At(iJet); iB2=iJet; bTag2=jetB2->BTag; }

      }

      if (jet->TauTag!=0) {

	if ((iT1==-1)) { jetTau1 = (Jet*) branchJet->At(iJet); iT1=iJet; tauCat1=hadron;}
	else if ((iT1!=iJet) && (jet->PT > jetTau1->PT)) { jetTau2 = (Jet*) branchJet->At(iT1); iT2=iT1; tauCat2=hadron; jetTau1 = (Jet*) branchJet->At(iJet); iT1=iJet; tauCat1=hadron;}
	else if ((iT1!=iJet) && (iT2==-1)) { jetTau2 = (Jet*) branchJet->At(iJet); iT2=iJet; tauCat2=hadron;}
	else if ((iT1!=iJet) && (jet->PT > jetTau2->PT)) { jetTau2 = (Jet*) branchJet->At(iJet); iT2=iJet; tauCat2=hadron;}

      }

    } // end reco jet loop

    if ((iT1==-1) || (iT2==-1)) {

      // get muonic taus
      for (Int_t iMuon=0; iMuon<branchMuon->GetEntries(); iMuon++) { // reco muon loop
	mu = (Muon*) branchMuon->At(iMuon);

	if (iT1==-1) { muTau = (Muon*) branchMuon->At(iMuon); iT1=iMuon; tauCat1=muon; }
	else if ((tauCat1!=muon)&&(iT2==-1)) { muTau = (Muon*) branchMuon->At(iMuon); iT2=iMuon; tauCat2=muon; }
	else if ( mu->PT > muTau->PT ) { 
	  muTau = (Muon*) branchMuon->At(iMuon); 
	  if ((tauCat1==muon)) iT1=iMuon; else iT2=iMuon;
	}
      }
      
      // get electronic taus
      for (Int_t iEle=0; iEle<branchElectron->GetEntries(); iEle++) { // reco ele loop
	ele = (Electron*) branchElectron->At(iEle);

	if (iT1==-1) { eleTau = (Electron*) branchElectron->At(iEle); iT1=iEle; tauCat1=electron; }
	else if ((tauCat1!=electron)&&(iT2==-1)) { eleTau = (Electron*) branchElectron->At(iEle); iT2=iEle; tauCat2=electron; }
	else if ( ele->PT > eleTau->PT ) { 
	  eleTau = (Electron*) branchElectron->At(iEle); 
	  if ((tauCat1==electron)) iT1=iEle; else iT2=iEle;
	}
      }
    }

    // look for extra jets
    for (Int_t iJet=0; iJet<branchJet->GetEntries(); iJet++) { // reconstructed jet loop
      jet = (Jet*) branchJet->At(iJet);

      if ((iJet==iB1) || (iJet==iB2)) continue;
      if ( ((tauCat1==hadron)&&(iJet==iT1)) || ((tauCat2==hadron)&&(iJet==iT2)) ) continue;

      if (iE==-1) { extraJet = (Jet*) branchJet->At(iJet); iE=iJet; }
      else if (jet->PT > extraJet->PT) { extraJet = (Jet*) branchJet->At(iJet); iE=iJet; }

    }

    // fill 4-vector for leading b-jet
    LorentzVector vRecoB1(0,0,0,0);
    if (bTag1==0) sRecoB1 = &nothing;
    else {
      vRecoB1.SetPt(jetB1->PT);
      vRecoB1.SetEta(jetB1->Eta);
      vRecoB1.SetPhi(jetB1->Phi);
      vRecoB1.SetM(jetB1->Mass);
      sRecoB1 = &vRecoB1;
    }

    // fill 4-vector for second b-jet
    LorentzVector vRecoB2(0,0,0,0);
    if (bTag2==0) sRecoB2 = &nothing;
    else {
      vRecoB2.SetPt(jetB2->PT);
      vRecoB2.SetEta(jetB2->Eta);
      vRecoB2.SetPhi(jetB2->Phi);
      vRecoB2.SetM(jetB2->Mass);
      sRecoB2 = &vRecoB2;
    }

    // fill 4-vector for leading tau
    LorentzVector vRecoTau1(0,0,0,0);
    if ( tauCat1 == 0 ) sRecoTau1 = &nothing;
    else if ( tauCat1 == hadron ) {
      vRecoTau1.SetPt(jetTau1->PT);
      vRecoTau1.SetEta(jetTau1->Eta);
      vRecoTau1.SetPhi(jetTau1->Phi);
      vRecoTau1.SetM(jetTau1->Mass);
      sRecoTau1 = &vRecoTau1;
    }
    else if ( tauCat1 == muon ) {
      vRecoTau1.SetPt(muTau->PT);
      vRecoTau1.SetEta(muTau->Eta);
      vRecoTau1.SetPhi(muTau->Phi);
      vRecoTau1.SetM(MUON_MASS);
      sRecoTau1 = &vRecoTau1;
    }
    else if ( tauCat1 == electron ) {
      vRecoTau1.SetPt(eleTau->PT);
      vRecoTau1.SetEta(eleTau->Eta);
      vRecoTau1.SetPhi(eleTau->Phi);
      vRecoTau1.SetM(ELE_MASS);
      sRecoTau1 = &vRecoTau1;
    }

    // fill 4-vector for second tau
    LorentzVector vRecoTau2(0,0,0,0);
    if ( tauCat2 == 0 ) sRecoTau2 = &nothing;
    else if ( tauCat2 == hadron ) {
      vRecoTau2.SetPt(jetTau2->PT);
      vRecoTau2.SetEta(jetTau2->Eta);
      vRecoTau2.SetPhi(jetTau2->Phi);
      vRecoTau2.SetM(jetTau2->Mass);
      sRecoTau2 = &vRecoTau2;
    }
    else if ( tauCat2 == muon ) {
      vRecoTau2.SetPt(muTau->PT);
      vRecoTau2.SetEta(muTau->Eta);
      vRecoTau2.SetPhi(muTau->Phi);
      vRecoTau2.SetM(MUON_MASS);
      sRecoTau2 = &vRecoTau2;
    }
    else if ( tauCat2 == electron ) {
      vRecoTau2.SetPt(eleTau->PT);
      vRecoTau2.SetEta(eleTau->Eta);
      vRecoTau2.SetPhi(eleTau->Phi);
      vRecoTau2.SetM(ELE_MASS);
      sRecoTau2 = &vRecoTau2;
    }

    // fill 4-vector for extra jet
    LorentzVector vRecoExtraJet(0,0,0,0);
    if (iE==-1) sRecoExtra = &nothing;
    else {
      vRecoExtraJet.SetPt(extraJet->PT);
      vRecoExtraJet.SetEta(extraJet->Eta);
      vRecoExtraJet.SetPhi(extraJet->Phi);
      vRecoExtraJet.SetM(extraJet->Mass);
      sRecoExtra = &vRecoExtraJet;
    }

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
      
    }

    LorentzVector vGenTau1(0,0,0,0);
    if ( iGenTau1 == -1) sGenTau1 = &nothing;
    else {
      vGenTau1.SetPt(genTau1->PT);
      vGenTau1.SetEta(genTau1->Eta);
      vGenTau1.SetPhi(genTau1->Phi);
      vGenTau1.SetM(genTau1->Mass);
      sGenTau1 = &vGenTau1;
    }

    LorentzVector vGenTau2(0,0,0,0);
    if ( iGenTau2 == -1) sGenTau2 = &nothing;
    else {
      vGenTau2.SetPt(genTau2->PT);
      vGenTau2.SetEta(genTau2->Eta);
      vGenTau2.SetPhi(genTau2->Phi);
      vGenTau2.SetM(genTau2->Mass);
      sGenTau2 = &vGenTau2;
    }

    LorentzVector vGenB1(0,0,0,0);
    if ( iGenB1 == -1) sGenB1 = &nothing;
    else {
      vGenB1.SetPt(genB1->PT);
      vGenB1.SetEta(genB1->Eta);
      vGenB1.SetPhi(genB1->Phi);
      vGenB1.SetM(genB1->Mass);
      sGenB1 = &vGenB1;
    }

    LorentzVector vGenB2(0,0,0,0);
    if ( iGenB2 == -1) sGenB2 = &nothing;
    else {
      vGenB2.SetPt(genB2->PT);
      vGenB2.SetEta(genB2->Eta);
      vGenB2.SetPhi(genB2->Phi);
      vGenB2.SetM(genB2->Mass);
      sGenB2 = &vGenB2;
    }

    // match generator level jets to generator particles
    for (Int_t iJet=0; iJet<branchGenJet->GetEntries(); iJet++) { // generator level jet loop
      genJet = (Jet*) branchGenJet->At(iJet);

      if ( (iGenB1 !=-1) && (deltaR(genJet->Eta, vGenB1.Eta(), genJet->Phi, vGenB1.Phi()) < MAX_MATCH_DIST) ) {
	iGenJetB1=iJet;
	genJetB1 = (Jet*) branchGenJet->At(iGenJetB1);
      }

      else if ( (iGenB2 !=-1) && (deltaR(genJet->Eta, vGenB2.Eta(), genJet->Phi, vGenB2.Phi()) < MAX_MATCH_DIST) ) {
	iGenJetB2=iJet;
	genJetB2 = (Jet*) branchGenJet->At(iGenJetB2);
      }

      else if ( (iGenTau1 !=-1) && (deltaR(genJet->Eta, vGenTau1.Eta(), genJet->Phi, vGenTau1.Phi()) < MAX_MATCH_DIST) ) {
	iGenJetTau1=iJet;
	genJetTau1 = (Jet*) branchGenJet->At(iGenJetTau1);
      }

      else if ( (iGenTau2 !=-1) && (deltaR(genJet->Eta, vGenTau2.Eta(), genJet->Phi, vGenTau2.Phi()) < MAX_MATCH_DIST) ) {
	iGenJetTau1=iJet;
	genJetTau2 = (Jet*) branchGenJet->At(iGenJetTau2);
      }

    }

    LorentzVector vGenJetTau1(0,0,0,0);
    if ( iGenJetTau1 == -1) sGenJetTau1 = &nothing;
    else {
      vGenJetTau1.SetPt(genJetTau1->PT);
      vGenJetTau1.SetEta(genJetTau1->Eta);
      vGenJetTau1.SetPhi(genJetTau1->Phi);
      vGenJetTau1.SetM(genJetTau1->Mass);
      sGenJetTau1 = &vGenJetTau1;
    }

    LorentzVector vGenJetTau2(0,0,0,0);
    if ( iGenJetTau2 == -1) sGenJetTau2 = &nothing;
    else {
      vGenJetTau2.SetPt(genJetTau2->PT);
      vGenJetTau2.SetEta(genJetTau2->Eta);
      vGenJetTau2.SetPhi(genJetTau2->Phi);
      vGenJetTau2.SetM(genJetTau2->Mass);
      sGenJetTau2 = &vGenJetTau2;
    }

    LorentzVector vGenJetB1(0,0,0,0);
    if ( iGenJetB1 == -1) sGenJetB1 = &nothing;
    else {
      vGenJetB1.SetPt(genJetB1->PT);
      vGenJetB1.SetEta(genJetB1->Eta);
      vGenJetB1.SetPhi(genJetB1->Phi);
      vGenJetB1.SetM(genJetB1->Mass);
      sGenJetB1 = &vGenJetB1;
    }

    LorentzVector vGenJetB2(0,0,0,0);
    if ( iGenJetB2 == -1) sGenJetB2 = &nothing;
    else {
      vGenJetB2.SetPt(genJetB2->PT);
      vGenJetB2.SetEta(genJetB2->Eta);
      vGenJetB2.SetPhi(genJetB2->Phi);
      vGenJetB2.SetM(genJetB2->Mass);
      sGenJetB2 = &vGenJetB2;
    }

    // ********************
    // MT2 CALC
    // ********************

    missET = (MissingET*) branchMET->At(0);

    met=missET->MET;
    metPhi=missET->Phi;

    if ( (iB1!=-1) && (iB2!=-1) && (iT1!=-1) && (iT2!=-1) ) {

      tau1.SetMagPhi(sRecoTau1->Pt(), sRecoTau1->Phi());
      tau2.SetMagPhi(sRecoTau2->Pt(), sRecoTau2->Phi());
      mTau1=sRecoTau1->M();
      mTau2=sRecoTau2->M();
      
      b1.SetMagPhi(sRecoB1->Pt(), sRecoB1->Phi());
      b2.SetMagPhi(sRecoB2->Pt(), sRecoB2->Phi());
      mB1=sRecoB1->M();
      mB2=sRecoB2->M();
      
      mpt.SetMagPhi(met, metPhi);
      
      TVector2 sumPt = tau1+tau2+mpt;
      
      smT2 testing = smT2();
      testing.SetB1(b1);
      testing.SetB2(b2);
      testing.SetMPT(sumPt);
      testing.SetMB1(mB1);
      testing.SetMB2(mB2);
      testing.SetMT1(mTau1);
      testing.SetMT2(mTau2);
      
      TVector2 c1=sumPt;
      TVector2 c2=sumPt-c1;
      
      ROOT::Math::Functor f(testing,2);
      double step[2] = {0.1, 0.1};
      double variable[2] = { 0.5*c1.Mod(), 0.0 };
      
      min->SetFunction(f);
      min->SetLimitedVariable(0,"cT",variable[0], step[0], 0.0, sumPt.Mod());
      min->SetLimitedVariable(1,"cPhi",variable[1], step[1], 0.0, TMath::Pi());
      
      min->Minimize();
      
      mt2 = min->MinValue();
    }
    
    else { mt2 = -1; }

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
