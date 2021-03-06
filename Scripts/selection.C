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

Int_t puJetID( Float_t eta, Float_t meanSqDeltaR, Float_t beta );

void selection(const TString inputfile="/afs/cern.ch/work/j/jlawhorn/public/HHToTTBB_14TeV/HHToTTBB_14TeV_0.root",
	       const Float_t xsec=2.92,
	       const Float_t totalEvents=5000,
	       Int_t   sampleNo=100,
	       const TString outputfile="test.root") {

  cout << inputfile << " " << xsec << " " << totalEvents << " " << outputfile << endl;

  // declare constants
  const Double_t MUON_MASS = 0.105658369;
  const Double_t ELE_MASS  = 0.000511;

  const Int_t TAU_ID_CODE = 15;
  const Int_t B_ID_CODE = 5;

  const Float_t MAX_MATCH_DIST = 0.5;

  // event types
  enum { HH=0, TT, H, EWK, ZJET, ETC };

  // tau decay modes
  enum { hadron=1, electron, muon };

  // setup mt2 minimizer
  ROOT::Math::Minimizer *min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
  min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
  min->SetTolerance(10.0);
  min->SetPrintLevel(0);

  TVector2 tau1(0,0), tau2(0,0), mpt(0,0), ppmpt(0,0);
  TVector2 b1(0,0), b2(0,0);
  Float_t mTau1=0, mTau2=0;
  Float_t mB1=0, mB2=0;
  Double_t mt2=0;
  Double_t ppMt2=0;

  // read input input file
  TChain chain("Delphes");
  chain.Add(inputfile);
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchMET =treeReader->UseBranch("MissingET");
  TClonesArray *branchPuppiMET =treeReader->UseBranch("PuppiMissingET");

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
  Jet *extraJet=0;

  Int_t iB1=-1,         iB2=-1;
  Int_t iT1=-1,         iT2=-1;
  Int_t iGenTau1=-1,    iGenTau2=-1;
  Int_t iGenB1=-1,      iGenB2=-1;
  Int_t iGenJetTau1=-1, iGenJetTau2=-1;
  Int_t iGenJetB1=-1,   iGenJetB2=-1;
  Int_t iExtra=-1;

  Int_t eventType;
  Int_t genInfo;
  Float_t eventWeight;
  Float_t met, metPhi;
  Float_t ppMet, ppMetPhi;
  Int_t tauCat1=0, tauCat2=0;
  Int_t bTag1=0, bTag2=0;

  Int_t tFake1=0, tFake2=0;
  Int_t bFake1=0, bFake2=0;

  Float_t ptTau1, ptTau2, ptB1, ptB2;
  Float_t etaTau1, etaTau2, etaB1, etaB2;
  Float_t phiTau1, phiTau2, phiB1, phiB2;

  Float_t mTT=0, mJJ=0, mHH=0, ptHH=0;

  LorentzVector *sRecoTau1=0, *sRecoTau2=0;
  LorentzVector *sGenJetTau1=0, *sGenJetTau2=0;
  LorentzVector *sGenTau1=0, *sGenTau2=0;

  LorentzVector *sRecoB1=0, *sRecoB2=0;
  LorentzVector *sGenJetB1=0, *sGenJetB2=0;
  LorentzVector *sGenB1=0, *sGenB2=0;

  LorentzVector *sRecoJet=0;

  TFile *outFile = new TFile(outputfile, "RECREATE");

  // tree to hold information about selected events
  TTree *outTree = new TTree("Events", "Events");
  outTree->Branch("eventWeight",    &eventWeight,    "eventWeight/f");  // event weight from cross-section and Event->Weight
  outTree->Branch("sampleNo",       &sampleNo,       "sampleNo/i");     // sample number (see config card for details)
  outTree->Branch("genInfo",        &genInfo,        "genInfo/i");      // generator level info (see below)
  outTree->Branch("eventType",      &eventType,      "eventType/i");    // see line 47
  outTree->Branch("tauCat1",        &tauCat1,        "tauCat1/i");      // leading tau final state - jet, muon, electron
  outTree->Branch("tauCat2",        &tauCat2,        "tauCat2/i");      // second tau final state - jet, muon, electron
  outTree->Branch("bTag1",          &bTag1,          "bTag1/i");        // leading b-jet tag from delphes
  outTree->Branch("bTag2",          &bTag2,          "bTag2/i");        // second b-jet tag from delphes
  outTree->Branch("tFake1",         &tFake1,         "tFake1/i");       // is tau1 matched to gen tau?
  outTree->Branch("tFake2",         &tFake2,         "tFake2/i");       // is tau2 matched to gen tau?
  outTree->Branch("bFake1",         &bFake1,         "bFake1/i");       // is b1 matched to gen b?
  outTree->Branch("bFake2",         &bFake2,         "bFake2/i");       // is b2 matched to gen b?
  outTree->Branch("ptTau1",         &ptTau1,         "ptTau1/f");       // pt(Tau1)
  outTree->Branch("etaTau1",        &etaTau1,        "etaTau1/f");      // eta(Tau1)
  outTree->Branch("phiTau1",        &phiTau1,        "phiTau1/f");      // phi(Tau1)
  outTree->Branch("ptTau2",         &ptTau2,         "ptTau2/f");       // pt(Tau2)
  outTree->Branch("etaTau2",        &etaTau2,        "etaTau2/f");      // eta(Tau2)
  outTree->Branch("phiTau2",        &phiTau2,        "phiTau2/f");      // phi(Tau2)
  outTree->Branch("ptB1",           &ptB1,           "ptB1/f");         // pt(B1)
  outTree->Branch("etaB1",          &etaB1,          "etaB1/f");        // eta(B1)
  outTree->Branch("phiB1",          &phiB1,          "phiB1/f");        // phi(B1)
  outTree->Branch("ptB2",           &ptB2,           "ptB2/f");         // pt(B2)
  outTree->Branch("etaB2",          &etaB2,          "etaB2/f");        // eta(B2)
  outTree->Branch("phiB2",          &phiB2,          "phiB2/f");        // phi(B2)
  outTree->Branch("met",            &met,            "met/f");          // missing transverse energy
  outTree->Branch("metPhi",         &metPhi,         "metPhi/f");       // missing transverse energy phi
  outTree->Branch("ppMet",          &ppMet,          "ppMet/f");        // PUPPI missing transverse energy
  outTree->Branch("ppMetPhi",       &ppMetPhi,       "ppMetPhi/f");     // PUPPI missing transverse energy phi
  outTree->Branch("mt2",            &mt2,            "mt2/D");          // "stransverse mass"
  outTree->Branch("ppMt2",          &ppMt2,          "ppMt2/D");        // PUPPI "stransverse mass"
  outTree->Branch("mTT",            &mTT,            "mTT/f");          // mass(tautau)
  outTree->Branch("mJJ",            &mJJ,            "mJJ/f");          // mass(jetjet)
  outTree->Branch("mHH",            &mHH,            "mHH/f");          // mass(HH)
  outTree->Branch("ptHH",           &ptHH,           "ptHH/f");         // pt(HH)
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
  outTree->Branch("sRecoJet",       "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sRecoJet);      // 4-vector for reconstructed extra jet

  // define placeholder vector for things that don't exist
  LorentzVector nothing(-999,-999,0,-999);

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
    iExtra=-1;
    met=0; metPhi=0; ppMet=0; ppMetPhi=0;
    tauCat1=0; tauCat2=0; bTag1=0; bTag2=0;
    tFake1=2; tFake2=2; bFake1=2; bFake2=2;
    mTT=-999; mJJ=-999; mHH=-999; ptHH=-999;
    ptTau1=-999; etaTau1=-999; phiTau1=-999;
    ptTau2=-999; etaTau2=-999; phiTau2=-999;
    ptB1=-999; etaB1=-999; phiB1=-999;
    ptB2=-999; etaB2=-999; phiB2=-999;
    eventType=-1;

    jetTau1=0; jetTau2=0; eleTau=0; muTau=0;
    jetB1=0;   jetB2=0; extraJet=0;
    genTau1=0; genTau2=0; genB1=0;  genB2=0;
    genJetB1=0; genJetB2=0; genJetTau1=0; genJetTau2=0;
    sGenTau1=0; sGenTau2=0; sGenB1=0;  sGenB2=0;
    sGenJetB1=0; sGenJetB2=0; sGenJetTau1=0; sGenJetTau2=0;
    sRecoJet=0;

    // ********************
    // EVENT WEIGHT
    // ********************

    // comment out following line for di-higgs samples
    //event = (LHEFEvent*) branchEvent->At(0);
    eventWeight = 1;
    if (branchEvent) {
      event = (LHEFEvent*) branchEvent->At(0);
      eventWeight*=event->Weight;
    }
    eventWeight *= xsec;
    eventWeight /= totalEvents;

    // ********************
    // RECO OBJECTS
    // ********************

    // get reconstructed hadronic taus
    for (Int_t iJet=0; iJet<branchJet->GetEntries(); iJet++) { // reconstructed jet loop
      jet = (Jet*) branchJet->At(iJet);

      if (fabs(jet->Eta)>4.0) continue;
      if (jet->PT<30) continue;

      if (puJetID(jet->Eta, jet->MeanSqDeltaR, jet->BetaStar)==1) continue;

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

      if (fabs(jet->Eta)>4.0) continue;
      if (jet->PT<30) continue;

      if (puJetID(jet->Eta, jet->MeanSqDeltaR, jet->BetaStar)==1) continue;

      if (jet->BTag==0) continue;
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
	jetB2 = (Jet*) branchJet->At(iB2); 
	bTag2=jetB2->BTag;
      }
    }

    for (Int_t iJet=0; iJet<branchJet->GetEntries(); iJet++) { // reconstructed jet loop
      jet = (Jet*) branchJet->At(iJet);

      if (fabs(jet->Eta)>4.0) continue;
      if (jet->PT<30) continue;

      if (puJetID(jet->Eta, jet->MeanSqDeltaR, jet->BetaStar)==1) continue;

      if ((jetTau1)&&(deltaR(jet->Eta, jetTau1->Eta, jet->Phi, jetTau1->Phi) < MAX_MATCH_DIST)) continue;
      if ((jetTau2)&&(deltaR(jet->Eta, jetTau2->Eta, jet->Phi, jetTau2->Phi) < MAX_MATCH_DIST)) continue;
      if ((jetB1)&&(deltaR(jet->Eta, jetB1->Eta, jet->Phi, jetB1->Phi) < MAX_MATCH_DIST)) continue;

      if (iB1==-1) {
	iB1=iJet; 
	jetB1 = (Jet*) branchJet->At(iB1); 
	bTag1=jetB1->BTag;
      }
      else if ((jet->PT > jetB1->PT) && ((bTag1==0)||(bTag1==-1))) {
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
	jetB2 = (Jet*) branchJet->At(iB2); 
	bTag2=jetB2->BTag;
      }
    }

    if ((iT1==-1) || (iT2==-1)) {
      // get muonic taus
      for (Int_t iMuon=0; iMuon<branchMuon->GetEntries(); iMuon++) { // reco muon loop
	mu = (Muon*) branchMuon->At(iMuon);

	if (fabs(mu->Eta)>4.0) continue;
	if (mu->PT<30) continue;

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

	if (fabs(ele->Eta)>4.0) continue;
	if (ele->PT<30) continue;

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

    for (Int_t iJet=0; iJet<branchJet->GetEntries(); iJet++) { // reconstructed jet loop
      jet = (Jet*) branchJet->At(iJet);

      if (fabs(jet->Eta)>4.0) continue;
      if (jet->PT<30) continue;

      if (puJetID(jet->Eta, jet->MeanSqDeltaR, jet->BetaStar)==1) continue;

      if ((jetTau1)&&(deltaR(jet->Eta, jetTau1->Eta, jet->Phi, jetTau1->Phi) < MAX_MATCH_DIST)) continue;
      if ((jetTau2)&&(deltaR(jet->Eta, jetTau2->Eta, jet->Phi, jetTau2->Phi) < MAX_MATCH_DIST)) continue;
      if ((jetB1)&&(deltaR(jet->Eta, jetB1->Eta, jet->Phi, jetB1->Phi) < MAX_MATCH_DIST)) continue;
      if ((jetB2)&&(deltaR(jet->Eta, jetB2->Eta, jet->Phi, jetB2->Phi) < MAX_MATCH_DIST)) continue;

      if (iExtra==-1) {
	iExtra=iJet;
	extraJet=(Jet*)branchJet->At(iExtra);
      }
      else if (jet->PT > extraJet->PT) {
	iExtra=iJet;
	extraJet=(Jet*)branchJet->At(iExtra);
      }
    }

    // fill 4-vector for leading b-jet
    LorentzVector vRecoB1(0,0,0,0);
    if (jetB1) {
      vRecoB1.SetPt(jetB1->PT);
      vRecoB1.SetEta(jetB1->Eta);
      vRecoB1.SetPhi(jetB1->Phi);
      vRecoB1.SetM(jetB1->Mass);
      sRecoB1 = &vRecoB1;
      ptB1=jetB1->PT;
      etaB1=jetB1->Eta;
      phiB1=jetB1->Phi;
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
      ptB2=jetB2->PT;
      etaB2=jetB2->Eta;
      phiB2=jetB2->Phi;
    }
    else sRecoB2 = &nothing;

    // fill 4-vector for leading tau
    LorentzVector vRecoTau1(0,0,0,0);
    if (jetTau1) {
      vRecoTau1.SetPt(jetTau1->PT);
      vRecoTau1.SetEta(jetTau1->Eta);
      vRecoTau1.SetPhi(jetTau1->Phi);
      vRecoTau1.SetM(jetTau1->Mass);
      sRecoTau1 = &vRecoTau1;
      ptTau1=jetTau1->PT;
      etaTau1=jetTau1->Eta;
      phiTau1=jetTau1->Phi;
    }
    else if ((muTau)&&(tauCat1==muon)) {
      vRecoTau1.SetPt(muTau->PT);
      vRecoTau1.SetEta(muTau->Eta);
      vRecoTau1.SetPhi(muTau->Phi);
      vRecoTau1.SetM(MUON_MASS);
      sRecoTau1 = &vRecoTau1;
      ptTau1=muTau->PT;
      etaTau1=muTau->Eta;
      phiTau1=muTau->Phi;
    }
    else if ((eleTau)&&(tauCat1==electron)) {
      vRecoTau1.SetPt(eleTau->PT);
      vRecoTau1.SetEta(eleTau->Eta);
      vRecoTau1.SetPhi(eleTau->Phi);
      vRecoTau1.SetM(ELE_MASS);
      sRecoTau1 = &vRecoTau1;
      ptTau1=eleTau->PT;
      etaTau1=eleTau->Eta;
      phiTau1=eleTau->Phi;
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
      ptTau2=jetTau2->PT;
      etaTau2=jetTau2->Eta;
      phiTau2=jetTau2->Phi;
    }
    else if ((muTau)&&(tauCat2==muon)) {
      vRecoTau2.SetPt(muTau->PT);
      vRecoTau2.SetEta(muTau->Eta);
      vRecoTau2.SetPhi(muTau->Phi);
      vRecoTau2.SetM(MUON_MASS);
      sRecoTau2 = &vRecoTau2;
      ptTau2=muTau->PT;
      etaTau2=muTau->Eta;
      phiTau2=muTau->Phi;
    }
    else if ((eleTau)&&(tauCat2==electron)) {
      vRecoTau2.SetPt(eleTau->PT);
      vRecoTau2.SetEta(eleTau->Eta);
      vRecoTau2.SetPhi(eleTau->Phi);
      vRecoTau2.SetM(ELE_MASS);
      sRecoTau2 = &vRecoTau2;
      ptTau2=eleTau->PT;
      etaTau2=eleTau->Eta;
      phiTau2=eleTau->Phi;
    }
    else sRecoTau2 = &nothing;

    LorentzVector vRecoJet(0,0,0,0);
    if (extraJet) {
      vRecoJet.SetPt(extraJet->PT);
      vRecoJet.SetEta(extraJet->Eta);
      vRecoJet.SetPhi(extraJet->Phi);
      vRecoJet.SetM(extraJet->Mass);
      sRecoJet = &vRecoJet;
    }
    else sRecoJet = &nothing;

    LorentzVector vTT = vRecoTau1+vRecoTau2;
    mTT=vTT.M();
    LorentzVector vJJ = vRecoB1+vRecoB2;
    mJJ=vJJ.M();
    LorentzVector vHH = vTT+vJJ;
    mHH=vHH.M();
    ptHH=vHH.Pt();

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
	else if ( (tauCat2 == -1) && (iGenTau1==-1) ) { 
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

    if ((sGenTau1->Pt()!=999)&&(tauCat1!=0)) { tFake1=0; }
    else if (tauCat1==1) { tFake1=1; }
    else {tFake1=2;}

    if ((sGenTau2->Pt()!=999)&&(tauCat2!=0)) { tFake2=0; }
    else if (tauCat2==1) { tFake2=1; }
    else {tFake2=2;}

    if ((sGenB1->Pt()!=999)&&(bTag1!=0)) { bFake1=0; }
    else if (bTag1!=0) { bFake1=1; }
    else {bFake1=2;}

    if ((sGenB2->Pt()!=999)&&(bTag2!=0)) { bFake2=0; }
    else if (bTag2!=0) { bFake2=1; }
    else {bFake2=2;}

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

    missET = (MissingET*) branchPuppiMET->At(0);
    
    ppMet=missET->MET;
    ppMetPhi=missET->Phi;

    if ( (sRecoTau1) && (sRecoTau2) && (sRecoB1) && (sRecoB2) ) {
    //if (0) {

      tau1.SetMagPhi(sRecoTau1->Pt(), sRecoTau1->Phi());
      tau2.SetMagPhi(sRecoTau2->Pt(), sRecoTau2->Phi());
      mTau1=sRecoTau1->M();
      mTau2=sRecoTau2->M();
      
      b1.SetMagPhi(sRecoB1->Pt(), sRecoB1->Phi());
      b2.SetMagPhi(sRecoB2->Pt(), sRecoB2->Phi());
      mB1=sRecoB1->M();
      mB2=sRecoB2->M();
      
      mpt.SetMagPhi(met, metPhi);
      ppmpt.SetMagPhi(ppMet, ppMetPhi);
      
      TVector2 sumPt = tau1+tau2+mpt;
      TVector2 ppSumPt = tau1+tau2+ppmpt;
      
      smT2 calcmt2 = smT2();
      calcmt2.SetB1(b1);
      calcmt2.SetB2(b2);
      calcmt2.SetMPT(sumPt);
      calcmt2.SetMB1(mB1);
      calcmt2.SetMB2(mB2);
      calcmt2.SetMT1(mTau1);
      calcmt2.SetMT2(mTau2);

      TVector2 c1=sumPt;
      TVector2 c2=sumPt-c1;
      
      ROOT::Math::Functor f(calcmt2,2);
      double step[2] = {0.1, 0.1};
      double variable[2] = { 0.5*c1.Mod(), 0.0 };
      
      min->SetFunction(f);
      min->SetLimitedVariable(0,"cT",variable[0], step[0], 0.0, sumPt.Mod());
      min->SetLimitedVariable(1,"cPhi",variable[1], step[1], 0.0, TMath::Pi());
      
      min->Minimize();
      mt2 = min->MinValue();

      calcmt2.SetMPT(ppSumPt);
      c1=ppSumPt;
      c2=sumPt-c1;

      ROOT::Math::Functor f2(calcmt2,2);
      step[0] = 0.1; step[1] = 0.1;
      variable[0] = 0.5*c1.Mod(); variable[1] = 0.0;

      min->SetFunction(f2);
      min->SetLimitedVariable(0,"cT",variable[0], step[0], 0.0, ppSumPt.Mod());
      min->SetLimitedVariable(1,"cPhi",variable[1], step[1], 0.0, TMath::Pi());
      
      min->Minimize();
      ppMt2 = min->MinValue();
            
    }
    
    else { 
      mt2 = 99999;
      ppMt2 = 99999;
    }

    // ********************
    // BKGD SORTING
    // ********************

    Int_t nH=0, nW=0, nZ=0, nT=0, nB=0, nQ=0, nG=0, nP=0, nL=0;
    for (Int_t iPart=0; iPart<branchParticle->GetEntries(); iPart++) {
      genParticle = (GenParticle*) branchParticle->At(iPart);
      
      Int_t pid=fabs(genParticle->PID);
      if (pid==2212) continue;
      if (pid==25) nH++;
      else if (pid==23) nZ++;
      else if (pid==24) nW++;
      else if (pid==6) nT++;
      else if (pid==5) nB++;
      else if (pid<5)  nQ++;
      else if (pid==21) nG++;
      else if (pid==22) nP++;
      else if (pid<17 && pid>10) nL++;
    }

    if (nH>10) nH=9; if (nW>10) nW=9;
    if (nZ>10) nZ=9; if (nT>10) nT=9;
    if (nB>10) nB=9; if (nQ>10) nQ=9;
    if (nG>10) nG=9; if (nP>10) nP=9;
    if (nL>10) nL=9;

    genInfo=nH+1e1*nW+1e2*nZ+1e3*nT+1e4*nB+1e5*nQ+1e6*nG+1e7*nP+1e8*nL;

    //cout << setprecision(9) << genInfo << endl;

    if ( nH==2 ) eventType=HH;
    else if ( nH>0 ) eventType=H;
    else if ( nT==2 ) eventType=TT;
    else if ( nZ==1 && (nG+nQ)>0 ) eventType=ZJET;
    else if ( nW+nZ+nT>0 ) eventType=EWK;
    else eventType=ETC;

    //cout << eventType << endl;

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

Int_t puJetID( Float_t eta, Float_t meanSqDeltaR, Float_t betastar) {
  
  Float_t MeanSqDeltaRMaxBarrel=0.07;
  Float_t BetaMinBarrel=0.87;
  Float_t MeanSqDeltaRMaxEndcap=0.07;
  Float_t BetaMinEndcap=0.85;

  //cout << eta << ", " << meanSqDeltaR << ", " << betastar << ": ";

  if (fabs(eta)<1.5) {
    if ((meanSqDeltaR<MeanSqDeltaRMaxBarrel)&&(betastar<BetaMinBarrel)) {
      //cout << "barrel 0" << endl;
      return 0;
    }
    else {
      //cout << "barrel 1" << endl;
      return 1;
    }
  }
  else if (fabs(eta)<4.0) {
    if ((meanSqDeltaR<MeanSqDeltaRMaxEndcap)&&(betastar<BetaMinEndcap)) {
      //cout << "endcap 0" << endl;
      return 0;
    }
    else {
      //cout << "endcap 1" << endl;
      return 1;
    }
  }
  //cout << "forward 1" << endl;
  return 1;

}
