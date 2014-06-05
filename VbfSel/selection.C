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

Float_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 );

Int_t puJetID( Float_t eta, Float_t meanSqDeltaR, Float_t beta );

void selection(const TString inputfile="root://eoscms.cern.ch//eos/cms/store/group/upgrade/delphes/test4/Bjj-vbf-4p-0-700-v1510_14TEV_152178247_PhaseII_Conf4_140PileUp.root",
	 const Float_t xsec=1.0,
	 const Float_t totalEvents=100,
	 const TString outputfile="test.root") {

  cout << inputfile << " " << xsec << " " << totalEvents << " " << outputfile << endl;

  // declare constants
  const Double_t MUON_MASS = 0.105658369;
  const Double_t ELE_MASS  = 0.000511;

  const Int_t TAU_ID_CODE = 15;

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
  Jet *jet1=0, *jet2=0;
  GenParticle *genTau1=0, *genTau2=0;
  Jet *genJetTau1=0, *genJetTau2=0;

  Int_t iJet1=-1,       iJet2=-1;
  Int_t iT1=-1,         iT2=-1;
  Int_t iGenTau1=-1,    iGenTau2=-1;
  Int_t iGenJetTau1=-1, iGenJetTau2=-1;

  Int_t nCentral=0;

  Int_t eventType;
  Float_t eventWeight;

  Float_t ptTau1, ptTau2, ptJet1, ptJet2;
  Float_t etaTau1, etaTau2, etaJet1, etaJet2;
  Float_t phiTau1, phiTau2, phiJet1, phiJet2;

  Float_t met, metPhi;
  Float_t dEta, mJJ, mTT;
  Int_t tauCat1=0, tauCat2=0;
  Int_t tFake1=0, tFake2=0;

  LorentzVector *sRecoTau1=0, *sRecoTau2=0;
  LorentzVector *sGenJetTau1=0, *sGenJetTau2=0;
  LorentzVector *sGenTau1=0, *sGenTau2=0;
  LorentzVector *sRecoJet1=0, *sRecoJet2=0;

  TFile *outFile = new TFile(outputfile, "RECREATE");

  // tree to hold information about selected events
  TTree *outTree = new TTree("Events", "Events");
  outTree->Branch("eventWeight",    &eventWeight,    "eventWeight/f");  // event weight from cross-section and Event->Weight
  outTree->Branch("eventType",      &eventType,      "eventType/i");    // event type (0=signal, 1=tt, 2=zh, 3=other)
  outTree->Branch("tauCat1",        &tauCat1,        "tauCat1/i");      // leading tau final state - jet, muon, electron
  outTree->Branch("tauCat2",        &tauCat2,        "tauCat2/i");      // second tau final state - jet, muon, electron
  outTree->Branch("tFake1",         &tFake1,         "tFake1/i");    
  outTree->Branch("tFake2",         &tFake2,         "tFake2/i");    
  outTree->Branch("ptTau1",         &ptTau1,         "ptTau1/f");       // pt(Tau1)                                                             
  outTree->Branch("etaTau1",        &etaTau1,        "etaTau1/f");      // eta(Tau1)                                                            
  outTree->Branch("phiTau1",        &phiTau1,        "phiTau1/f");      // phi(Tau1)                                                            
  outTree->Branch("ptTau2",         &ptTau2,         "ptTau2/f");       // pt(Tau2)                                                             
  outTree->Branch("etaTau2",        &etaTau2,        "etaTau2/f");      // eta(Tau2)                                                            
  outTree->Branch("phiTau2",        &phiTau2,        "phiTau2/f");      // phi(Tau2)                                                            
  outTree->Branch("ptJet1",         &ptJet1,         "ptJet1/f");       // pt(Jet1)                                                               
  outTree->Branch("etaJet1",        &etaJet1,        "etaJet1/f");      // eta(Jet1)                                                              
  outTree->Branch("phiJet1",        &phiJet1,        "phiJet1/f");      // phi(Jet1)                                                              
  outTree->Branch("ptJet2",         &ptJet2,         "ptJet2/f");       // pt(Jet2)                                                               
  outTree->Branch("etaJet2",        &etaJet2,        "etaJet2/f");      // eta(Jet2)                                                              
  outTree->Branch("phiJet2",        &phiJet2,        "phiJet2/f");      // phi(Jet2)          
  outTree->Branch("met",            &met,            "met/f");          // missing transverse energy
  outTree->Branch("metPhi",         &metPhi,         "metPhi/f");       // missing transverse energy phi
  outTree->Branch("dEta",           &dEta,           "dEta/f");         // delta Eta
  outTree->Branch("mJJ",            &mJJ,            "mJJ/f") ;         // mass(jets)
  outTree->Branch("mTT",            &mTT,            "mTT/f") ;         // mass(tautau)
  outTree->Branch("sGenTau1",       "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sGenTau1);      // 4-vector for generator leading tau
  outTree->Branch("sGenTau2",       "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sGenTau2);      // 4-vector for generator second tau
  outTree->Branch("sGenJetTau1",    "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sGenJetTau1);   // 4-vector for generator leading tau
  outTree->Branch("sGenJetTau2",    "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sGenJetTau2);   // 4-vector for generator second tau
  outTree->Branch("sRecoTau1",      "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sRecoTau1);     // 4-vector for reconstructed leading tau
  outTree->Branch("sRecoTau2",      "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sRecoTau2);     // 4-vector for reconstructed second tau
  outTree->Branch("sRecoJet1",      "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sRecoJet1);     // 4-vector for reconstructed leading forward-jet
  outTree->Branch("sRecoJet2",      "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sRecoJet2);     // 4-vector for reconstructed second forward-jet

  // define placeholder vector for things that don't exist
  LorentzVector nothing(-999,-999,0,-999);

  for (Int_t iEntry=0; iEntry<numberOfEntries; iEntry++) { // entry loop
    treeReader->ReadEntry(iEntry);

    // ********************
    // RESET
    // ********************
    
    iJet1=-1; iJet2=-1; iT1=-1; iT2=-1;
    iGenTau1=-1;    iGenTau2=-1;
    iGenJetTau1=-1; iGenJetTau2=-1;
    tauCat1=-1; tauCat2=-1;
    tFake1=2; tFake2=2;
    eventType=-1;

    mTT=-999; mJJ=-999; dEta=-999;
    ptTau1=-999; etaTau1=-999; phiTau1=-999;
    ptTau2=-999; etaTau2=-999; phiTau2=-999;
    ptJet1=-999; etaJet1=-999; phiJet1=-999;
    ptJet2=-999; etaJet2=-999; phiJet2=-999;

    jetTau1=0; jetTau2=0; eleTau=0; muTau=0;
    jet1=0;   jet2=0;
    genTau1=0; genTau2=0;
    genJetTau1=0; genJetTau2=0;
    sGenTau1=0; sGenTau2=0;
    sRecoJet1=0; sRecoJet2=0;

    // ********************
    // EVENT WEIGHT
    // ********************

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

      if (jet->TauTag==0) continue;
      if (fabs(jet->Eta)>4.0) continue;
      if (jet->PT<30) continue;

      if (puJetID(jet->Eta, jet->MeanSqDeltaR, jet->BetaStar)==1) continue;

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

    // get VBF jets
    for (Int_t iJet=0; iJet<branchJet->GetEntries(); iJet++) { // reconstructed jet loop
      jet = (Jet*) branchJet->At(iJet);

      if (fabs(jet->Eta)>4.7) continue;
      if (jet->PT<30) continue;

      if (puJetID(jet->Eta, jet->MeanSqDeltaR, jet->BetaStar)==1) continue;

      if ((jet1)&&(deltaR(jet->Eta, jet1->Eta, jet->Phi, jet1->Phi) < MAX_MATCH_DIST)) continue;
      if ((jet2)&&(deltaR(jet->Eta, jet2->Eta, jet->Phi, jet2->Phi) < MAX_MATCH_DIST)) continue;
      if ((jetTau1)&&(deltaR(jet->Eta, jetTau1->Eta, jet->Phi, jetTau1->Phi) < MAX_MATCH_DIST)) continue;
      if ((jetTau2)&&(deltaR(jet->Eta, jetTau2->Eta, jet->Phi, jetTau2->Phi) < MAX_MATCH_DIST)) continue;
      
      if (iJet1==-1) {
	iJet1=iJet; 
	jet1 = (Jet*) branchJet->At(iJet1); 
      }
      else if (jet->PT > jet1->PT) {
	iJet2=iJet1; 
	jet2 = (Jet*) branchJet->At(iJet2); 
	iJet1=iJet;
	jet1 = (Jet*) branchJet->At(iJet1);
      }
      else if (iJet2==-1) { 
	iJet2=iJet; 
	jet2 = (Jet*) branchJet->At(iJet2); 
      }
      else if (jet->PT > jet2->PT) { 
	iJet2=iJet; 
	jet2 = (Jet*) branchJet->At(iJet2); 
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

    if ((iT1==-1)||(iT2==-1)||(iJet1==-1)||(iJet2==-1)) continue;

    nCentral=0;

    /*    cout << endl;
    cout << "New event" << endl;
    cout << endl;
    if (jetTau1 && jetTau2) cout << jetTau1->Eta << " " << jetTau1->PT << ", " << jetTau2->Eta << " " << jetTau2->PT << endl;
    cout << jet1->Eta << " " << jet1->PT << ", " << jet2->Eta << " " << jet2->PT << endl;
    cout << "---" << endl;*/

    for (Int_t iJet=0; iJet<branchJet->GetEntries(); iJet++) { // reconstructed jet loop
      jet = (Jet*) branchJet->At(iJet);

      if ((jetTau1)&&(deltaR(jet->Eta, jetTau1->Eta, jet->Phi, jetTau1->Phi) < MAX_MATCH_DIST)) continue;
      if ((jetTau2)&&(deltaR(jet->Eta, jetTau2->Eta, jet->Phi, jetTau2->Phi) < MAX_MATCH_DIST)) continue;

      if (fabs(jet->Eta)>4.7) continue;
      if (jet->PT<30) continue;

      if ( (jet1->Eta > jet2->Eta) && (jet->Eta > jet2->Eta) && (jet1->Eta > jet->Eta) ) {
	nCentral++;
	//cout << "central! " << jet->PT << " " << jet->Eta << endl;
      }
      else if ( (jet2->Eta > jet1->Eta) && (jet->Eta > jet1->Eta) && (jet2->Eta > jet->Eta) ) {
	nCentral++;
	//cout << "central! " << jet->PT << " " << jet->Eta << endl;
      }
    }
      
    if (nCentral>0) continue;

    // fill 4-vector for leading jet
    LorentzVector vReco1(0,0,0,0);
    if (jet1) {
      vReco1.SetPt(jet1->PT);
      vReco1.SetEta(jet1->Eta);
      vReco1.SetPhi(jet1->Phi);
      vReco1.SetM(jet1->Mass);
      sRecoJet1 = &vReco1;
      ptJet1=jet1->PT;
      etaJet1=jet1->Eta;
      phiJet1=jet1->Phi;
    }
    else sRecoJet1 = &nothing;
    // fill 4-vector for second jet
    LorentzVector vReco2(0,0,0,0);
    if (jet2) {
      vReco2.SetPt(jet2->PT);
      vReco2.SetEta(jet2->Eta);
      vReco2.SetPhi(jet2->Phi);
      vReco2.SetM(jet2->Mass);
      sRecoJet2 = &vReco2;
      ptJet2=jet2->PT;
      etaJet2=jet2->Eta;
      phiJet2=jet2->Phi;
    }
    else sRecoJet2 = &nothing;

    if (!jet1 || !jet2) continue;
    LorentzVector vJJ = vReco1+vReco2;
    mJJ=vJJ.M();

    dEta=jet2->Eta-jet1->Eta;

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

    LorentzVector vTT = vRecoTau1+vRecoTau2;
    mTT=vTT.M();

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

    if ((sGenTau1->Pt()!=999)&&(tauCat1!=0)) { tFake1=0; }
    else if (tauCat1==1) { tFake1=1; }
    else {tFake1=2;}

    if ((sGenTau2->Pt()!=999)&&(tauCat2!=0)) { tFake2=0; }
    else if (tauCat2==1) { tFake2=1; }
    else {tFake2=2;}

    // match generator level jets to generator particles
    for (Int_t iJet=0; iJet<branchGenJet->GetEntries(); iJet++) { // generator level jet loop
      genJet = (Jet*) branchGenJet->At(iJet);

      if ((genTau1) && (deltaR(genJet->Eta, genTau1->Eta, genJet->Phi, genTau1->Phi) < MAX_MATCH_DIST) ) {
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
