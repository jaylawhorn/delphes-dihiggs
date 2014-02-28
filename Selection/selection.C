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

  // event types
  enum { HH=0, TT, ZH, WH, WW, ZZ, ZW, ETC };

  const Float_t MAX_MATCH_DIST = 0.3;

  // tau decay modes
  enum { hadron=1, electron, muon };

  // setup mt2 minimizer
  ROOT::Math::Minimizer *min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
  // set tolerance , etc...
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

  //TClonesArray *branchGenJet = treeReader->UseBranch("GenJet");
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  //TClonesArray *branchEvent = treeReader->UseBranch("Event");

  //set up loop variables
  GenParticle *genParticle;
  //Jet *genJet;
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
  Jet *leadingJet;

  Int_t iLead=-1, iB1=-1, iB2=-1, iT1=-1, iT2=-1, iE=-1;

  Int_t nExtra=0, nSelect=0;

  // set up output variables and file
  Int_t nEvents;
  Int_t eventType;
  Float_t eventWeight;
  Float_t met, metPhi;
  UInt_t tauCat1=0, tauCat2=0;
  UInt_t bTag1=0, bTag2=0;
  LorentzVector *recoTau1=0, *recoTau2=0, *recoB1=0, *recoB2=0, *recoLeadJet=0, *recoExtraJet=0;
  //LorentzVector *genTau1=0, *genTau2=0, *genDecayTau1=0, *genDecayTau2=0, *recoTau1=0, *recoTau2=0;
  //LorentzVector *genB1=0, *genB2=0, *recoB1=0, *recoB2=0, *boostJet=0, *genBoostJet=0;

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
  outTree->Branch("recoTau1",       "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &recoTau1);     // 4-vector for reconstructed leading tau
  outTree->Branch("recoTau2",       "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &recoTau2);     // 4-vector for reconstructed second tau
  outTree->Branch("recoB1",         "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &recoB1);       // 4-vector for reconstructed leading b-jet
  outTree->Branch("recoB2",         "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &recoB2);       // 4-vector for reconstructed second b-jet
  outTree->Branch("recoLeadJet",    "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &recoLeadJet);  // 4-vector for reconstructed leading jet
  outTree->Branch("recoExtraJet",   "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &recoExtraJet); // 4-vector for reconstructed extra jet

  // define placeholder vector for things that don't exist
  LorentzVector nothing(999,999,0,999);

  for (Int_t iEntry=0; iEntry<numberOfEntries; iEntry++) { // entry loop
    treeReader->ReadEntry(iEntry);
    //cout << endl << "new event" << endl;
    iLead=-1; iB1=-1; iB2=-1; iT1=-1; iT2=-1; iE=-1;
    tauCat1=0; tauCat2=0; bTag1=0; bTag2=0;
    eventType=-1;

    // comment out following line for di-higgs samples
    //event = (LHEFEvent*) branchEvent->At(0);
    eventWeight = 1;
    eventWeight *= xsec;
    // comment out following line for di-higgs samples
    //eventWeight *= event->Weight;

    // get reconstructed hadronic taus and b jets
    for (Int_t iJet=0; iJet<branchJet->GetEntries(); iJet++) { // reconstructed jet loop
      jet = (Jet*) branchJet->At(iJet);

      //if (jet->PT < 30) continue;
      //if (fabs(jet->Eta) > 2.5) continue;

      if (iLead==-1) {leadingJet = (Jet*) branchJet->At(iJet); iLead=iJet; }
      else if (jet->PT > leadingJet->PT) { leadingJet = (Jet*) branchJet->At(iJet); iLead=iJet; }

      if (jet->BTag!=0) {
	if ((iB1==-1)) { jetB1 = (Jet*) branchJet->At(iJet); iB1=iJet; bTag1=jetB1->BTag;}
	else if ((iB1!=iJet) && (jet->PT > jetB1->PT)) { jetB2 = (Jet*) branchJet->At(iB1); iB2=iB1; bTag2=bTag1; jetB1 = (Jet*) branchJet->At(iJet); iB1=iJet; bTag1=jetB1->BTag;}
	else if ((iB1!=iJet) && (iB2==-1)) { jetB2 = (Jet*) branchJet->At(iJet); iB2=iJet; bTag2=jetB2->BTag;}
	else if ((iB1!=iJet) && (jet->PT > jetB2->PT)) { jetB2 = (Jet*) branchJet->At(iJet); iB2=iJet; bTag2=jetB2->BTag;}
      }

      if ((iJet==iB1) || (iJet==iB2)) continue;

      if (jet->TauTag!=0) {
	if ((iT1==-1)) { jetTau1 = (Jet*) branchJet->At(iJet); iT1=iJet; tauCat1=hadron;}
	else if ((iT1!=iJet) && (jet->PT > jetTau1->PT)) { jetTau2 = (Jet*) branchJet->At(iT1); iT2=iT1; tauCat2=hadron; jetTau1 = (Jet*) branchJet->At(iJet); iT1=iJet; tauCat1=hadron;}
	else if ((iT1!=iJet) && (iT2==-1)) { jetTau2 = (Jet*) branchJet->At(iJet); iT2=iJet; tauCat2=hadron;}
	else if ((iT1!=iJet) && (jet->PT > jetTau2->PT)) { jetTau2 = (Jet*) branchJet->At(iJet); iT2=iJet; tauCat2=hadron;}
      }

    } // end reco jet loop

    //if ((iB1==-1) || (iB2==-1)) continue;

    if ((iT1==-1) || (iT2==-1)) {

      // get muonic taus
      for (Int_t iMuon=0; iMuon<branchMuon->GetEntries(); iMuon++) { // reco muon loop
	mu = (Muon*) branchMuon->At(iMuon);
	//if (mu->PT < 20) continue;
	//if (fabs(mu->Eta) > 2.5) continue;	

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
	//if (ele->PT < 20) continue;
	//if (fabs(ele->Eta) > 2.5) continue;

	if (iT1==-1) { eleTau = (Electron*) branchElectron->At(iEle); iT1=iEle; tauCat1=electron; }
	else if ((tauCat1!=electron)&&(iT2==-1)) { eleTau = (Electron*) branchElectron->At(iEle); iT2=iEle; tauCat2=electron; }
	else if ( ele->PT > eleTau->PT ) { 
	  eleTau = (Electron*) branchElectron->At(iEle); 
	  if ((tauCat1==electron)) iT1=iEle; else iT2=iEle;
	}
      }
    }

    if ( (iT1==-1) || (iT2==-1) ) continue;

    // look for extra jets
    for (Int_t iJet=0; iJet<branchJet->GetEntries(); iJet++) { // reconstructed jet loop
      jet = (Jet*) branchJet->At(iJet);

      //if (jet->PT < 30) continue;
      //if (fabs(jet->Eta) > 2.5) continue;

      if ((iJet==iB1) || (iJet==iB2)) continue;
      if ( ((tauCat1==hadron)&&(iJet==iT1)) || ((tauCat2==hadron)&&(iJet==iT2)) ) continue;

      if (iE==-1) { extraJet = (Jet*) branchJet->At(iJet); iE=iJet; }
      else if (jet->PT > extraJet->PT) { extraJet = (Jet*) branchJet->At(iJet); iE=iJet; }

    }

    // fill 4-vector for leading b-jet
    LorentzVector vRecoB1(0,0,0,0);
    if (bTag1==0) recoB1 = &nothing;
    else {
      vRecoB1.SetPt(jetB1->PT);
      vRecoB1.SetEta(jetB1->Eta);
      vRecoB1.SetPhi(jetB1->Phi);
      vRecoB1.SetM(jetB1->Mass);
      recoB1 = &vRecoB1;
    }

    // fill 4-vector for second b-jet
    LorentzVector vRecoB2(0,0,0,0);
    if (bTag2==0) recoB2 = &nothing;
    else {
      vRecoB2.SetPt(jetB2->PT);
      vRecoB2.SetEta(jetB2->Eta);
      vRecoB2.SetPhi(jetB2->Phi);
      vRecoB2.SetM(jetB2->Mass);
      recoB2 = &vRecoB2;
    }

    // fill 4-vector for leading tau
    LorentzVector vRecoTau1(0,0,0,0);
    if ( tauCat1 == 0 ) recoTau1 = &nothing;
    else if ( tauCat1 == hadron ) {
      vRecoTau1.SetPt(jetTau1->PT);
      vRecoTau1.SetEta(jetTau1->Eta);
      vRecoTau1.SetPhi(jetTau1->Phi);
      vRecoTau1.SetM(jetTau1->Mass);
      recoTau1 = &vRecoTau1;
    }
    else if ( tauCat1 == muon ) {
      vRecoTau1.SetPt(muTau->PT);
      vRecoTau1.SetEta(muTau->Eta);
      vRecoTau1.SetPhi(muTau->Phi);
      vRecoTau1.SetM(MUON_MASS);
      recoTau1 = &vRecoTau1;
    }
    else if ( tauCat1 == electron ) {
      vRecoTau1.SetPt(eleTau->PT);
      vRecoTau1.SetEta(eleTau->Eta);
      vRecoTau1.SetPhi(eleTau->Phi);
      vRecoTau1.SetM(ELE_MASS);
      recoTau1 = &vRecoTau1;
    }

    // fill 4-vector for second tau
    LorentzVector vRecoTau2(0,0,0,0);
    if ( tauCat2 == 0 ) recoTau2 = &nothing;
    else if ( tauCat2 == hadron ) {
      vRecoTau2.SetPt(jetTau2->PT);
      vRecoTau2.SetEta(jetTau2->Eta);
      vRecoTau2.SetPhi(jetTau2->Phi);
      vRecoTau2.SetM(jetTau2->Mass);
      recoTau2 = &vRecoTau2;
    }
    else if ( tauCat2 == muon ) {
      vRecoTau2.SetPt(muTau->PT);
      vRecoTau2.SetEta(muTau->Eta);
      vRecoTau2.SetPhi(muTau->Phi);
      vRecoTau2.SetM(MUON_MASS);
      recoTau2 = &vRecoTau2;
    }
    else if ( tauCat2 == electron ) {
      vRecoTau2.SetPt(eleTau->PT);
      vRecoTau2.SetEta(eleTau->Eta);
      vRecoTau2.SetPhi(eleTau->Phi);
      vRecoTau2.SetM(ELE_MASS);
      recoTau2 = &vRecoTau2;
    }

    // fill 4-vector for leading jet
    LorentzVector vRecoLeadJet(leadingJet->PT,leadingJet->Eta,leadingJet->Phi,leadingJet->Mass);
    recoLeadJet = &vRecoLeadJet;

    nSelect++;

    // fill 4-vector for extra jet
    LorentzVector vRecoExtraJet(0,0,0,0);
    if (iE==-1) recoExtraJet = &nothing;
    else {
      nExtra++;
      vRecoExtraJet.SetPt(extraJet->PT);
      vRecoExtraJet.SetEta(extraJet->Eta);
      vRecoExtraJet.SetPhi(extraJet->Phi);
      vRecoExtraJet.SetM(extraJet->Mass);
      recoExtraJet = &vRecoExtraJet;
    }

    missET = (MissingET*) branchMET->At(0);

    met=missET->MET;
    metPhi=missET->Phi;

    tau1.SetMagPhi(recoTau1->Pt(), recoTau1->Phi());
    tau2.SetMagPhi(recoTau2->Pt(), recoTau2->Phi());
    mTau1=recoTau1->M();
    mTau2=recoTau2->M();

    b1.SetMagPhi(recoB1->Pt(), recoB1->Phi());
    b2.SetMagPhi(recoB2->Pt(), recoB2->Phi());
    mB1=recoB1->M();
    mB2=recoB2->M();

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

    //testing(c1.Mod(), c1.Phi());
    ROOT::Math::Functor f(testing,2);
    double step[2] = {0.1, 0.1};
    double variable[2] = { 0.5*c1.Mod(), 0.0 };

    min->SetFunction(f);
    min->SetLimitedVariable(0,"cT",variable[0], step[0], 0.0, sumPt.Mod());
    min->SetLimitedVariable(1,"cPhi",variable[1], step[1], 0.0, TMath::Pi());

    min->Minimize();

    mt2 = min->MinValue();
    //cout << "mt2: " << mt2 << endl;

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

    cout <<"H" << nH << ", W" << nW << ", Z" << nZ << ", T" << nT << " : " << eventType << endl;

    outTree->Fill();

  } // end event loop

  outFile->Write();
  outFile->Save();

  cout << nExtra << " events with an extra jet out of " << nSelect << " selected events" << endl;

}

Float_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 ) {

  const Float_t pi = 3.14159265358979;

  Float_t etaDiff = (eta1-eta2);
  Float_t phiDiff = fabs(phi1-phi2);
  while ( phiDiff>pi ) phiDiff = fabs(phiDiff-2.0*pi);

  Float_t deltaRSquared = etaDiff*etaDiff + phiDiff*phiDiff;

  return TMath::Sqrt(deltaRSquared);

}
