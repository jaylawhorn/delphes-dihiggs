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
#include <TProfile.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "Math/LorentzVector.h"
#include "TLorentzVector.h"
#include "modules/Delphes.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "classes/DelphesClasses.h"

#endif

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;

Float_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 );

//void selection(const TString infile="/afs/cern.ch/work/j/jlawhorn/public/delphes-test/test_B_all.root") {
//void selection(const TString infile="root://eoscms.cern.ch//store/group/upgrade/delphes/PhaseII_140PU_ProdJul28/B-4p-0-1-v1510_14TEV/B-4p-0-1-v1510_14TEV_92673021_PhaseII_Conf4_140PileUp_seed92673026_5of5.root") {
//void selection(const TString infile="root://eoscms.cern.ch//store/group/upgrade/delphes/alexis_tests/PhaseI_50PU/tt-4p-1700-2500-v1510_14TEV/tt-4p-1700-2500-v1510_14TEV_187877707_PhaseI_Conf0_50PileUp_seed187877711_4of5.root") {
void selection() { 
//const TString infile="root://eoscms.cern.ch//store/group/upgrade/delphes/PhaseII_140PU_ProdJul28/tt-4p-1700-2500-v1510_14TEV/tt-4p-1700-2500-v1510_14TEV_187877707_PhaseII_Conf4_140PileUp_seed187877711_4of5.root") {
//void selection(const TString infile="root://eoscms.cern.ch//store/group/phys_higgs/upgrade/PhaseII/Configuration4v2/140PileUp/tt-4p-0-600-v1510_14TEV/tt-4p-0-600-v1510_14TEV_96138438_PhaseII_Conf4v2_140PileUp.root") {
//void selection(const TString infile="root://eoscms.cern.ch//store/group/upgrade/delphes/ProdJun14/tt-4p-0-600-v1510_14TEV/tt-4p-0-600-v1510_14TEV_216084397_PhaseII_Conf4_140PileUp_seed216084401_4of5.root") {

  const Float_t MAX_MATCH_DIST = 0.4;

  TProfile *hEffPtEle = new TProfile("hEffPtEle", "hEffPtEle", 4, 0, 200);
  TProfile *hEffEtaEle = new TProfile("hEffEtaEle", "hEffEtaEle", 4, -4.0, 4.0);                                                                
  TProfile *hEffPtMu = new TProfile("hEffPtMu", "hEffPtMu", 4, 0, 200);                                                                    
  TProfile *hEffEtaMu = new TProfile("hEffEtaMu", "hEffEtaMu", 4, -4.0, 4.0);                                                                
  TProfile *hEffPtTau = new TProfile("hEffPtTau", "hEffPtTau", 4, 0, 200);                                                                    
  TProfile *hEffEtaTau = new TProfile("hEffEtaTau", "hEffEtaTau", 4, -4.0, 4.0);                                                                
  TProfile *hEffPtBJet = new TProfile("hEffPtBJet", "hEffPtBJet", 4, 0, 200);                                                                    
  TProfile *hEffEtaBJet = new TProfile("hEffEtaBJet", "hEffEtaBJet", 4, -4.0, 4.0);                                                                

  TProfile *hResPtEle = new TProfile("hResPtEle", "hResPtEle", 4, 0, 200);
  TProfile *hResEtaEle = new TProfile("hResEtaEle", "hResEtaEle", 4, -4.0, 4.0);
  TProfile *hResPtMu = new TProfile("hResPtMu", "hResPtMu", 4, 0, 200);
  TProfile *hResEtaMu = new TProfile("hResEtaMu", "hResEtaMu", 4, -4.0, 4.0);
  TProfile *hResPtTau = new TProfile("hResPtTau", "hResPtTau", 4, 0, 200);
  TProfile *hResEtaTau = new TProfile("hResEtaTau", "hResEtaTau", 4, -4.0, 4.0);
  TProfile *hResPtBJet = new TProfile("hResPtBJet", "hResPtBJet", 4, 0, 200);
  TProfile *hResEtaBJet = new TProfile("hResEtaBJet", "hResEtaBJet", 4, -4.0, 4.0);

  // tau decay modes
  enum { hadron=1, electron, muon };

  // read input input file
  TChain chain("Delphes");
  chain.Add("root://eoscms.cern.ch//store/group/upgrade/delphes/alexis_tests/PhaseI_50PU/tt-4p-1700-2500-v1510_14TEV/tt-4p-1700-2500-v1510_14TEV_186620505_PhaseI_Conf0_50PileUp_seed186620506_1of5.root");
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  TClonesArray *branchJet = treeReader->UseBranch("Jet");

  if (!(branchJet)) {
    cout << "File broken" << endl;
    return;
  }

  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
  TClonesArray *branchMET =treeReader->UseBranch("MissingET");
  TClonesArray *branchPuppiMET =treeReader->UseBranch("PuppiMissingET");
  TClonesArray *branchPileupMET =treeReader->UseBranch("PileUpJetIDMissingET");
  TClonesArray *branchRho = treeReader->UseBranch("Rho");
  TClonesArray *branchGenJet = treeReader->UseBranch("GenJet");
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchEvent = treeReader->UseBranch("Event");

  //set up loop variables
  GenParticle *genParticle=0;
  Jet *genJet=0;
  Jet *jet=0;
  Electron *ele=0;
  Muon *mu=0;
  Photon *gam=0;
  MissingET *missET=0;
  LHEFEvent *event=0;

  Int_t iEleMatch=-1, iMuMatch=-1, iTauMatch=-1, iBMatch=-1;

  Int_t nE=0, nM=0, nT=0, nB=0;
  Int_t mE=0, mM=0, mT=0, mB=0;
  Int_t mTauTag=0, mBTag=0;
  Int_t mTauTag1=0, mBTag1=0, mBTag2=0, mBTag3=0;
  Int_t nJ=0;

  for (Int_t iEntry=0; iEntry<numberOfEntries; iEntry++) { // entry loop
  //for (Int_t iEntry=0; iEntry<1000; iEntry++) { // entry loop
    treeReader->ReadEntry(iEntry);

    // ********************
    // GEN PARTICLES
    // ********************

    //cout << " ----- " << endl;

    for (Int_t iParticle=0; iParticle<branchParticle->GetEntries(); iParticle++) { // generator particle loop
      genParticle = (GenParticle*) branchParticle->At(iParticle);

      iEleMatch=-1; iMuMatch=-1; iTauMatch=-1; iBMatch=-1;
      
      if ( genParticle->Status!=3 || genParticle->PID==2212 ) continue;
      
      Int_t pid=fabs(genParticle->PID);
      
      if (pid==11) {
	nE++;
	for (Int_t iEle=0; iEle<branchElectron->GetEntries(); iEle++) { // electron loop
	  ele = (Electron*) branchElectron->At(iEle);
	  if ( deltaR(ele->Eta, genParticle->Eta, ele->Phi, genParticle->Phi) < MAX_MATCH_DIST ) {
	    mE++;
	    iEleMatch=iEle;
	  }
	}
	
	if (iEleMatch!=-1 && iEleMatch<branchElectron->GetEntries()) {
	  ele = (Electron*) branchElectron->At(iEleMatch);
	  hEffPtEle->Fill(genParticle->PT, 1);
	  hEffEtaEle->Fill(genParticle->Eta, 1);
	  hResPtEle->Fill(genParticle->PT, (ele->PT-genParticle->PT)/genParticle->PT);
	  hResEtaEle->Fill(genParticle->Eta, (ele->PT-genParticle->PT)/genParticle->PT);
	}
	else {
	  hEffPtEle->Fill(genParticle->PT, 0);
	  hEffEtaEle->Fill(genParticle->Eta, 0);
	}
      }
      if (pid==13) {
	nM++;
	for (Int_t iMu=0; iMu<branchMuon->GetEntries(); iMu++) { // muon loop
	  mu = (Muon*) branchMuon->At(iMu);
	  if ( deltaR(mu->Eta, genParticle->Eta, mu->Phi, genParticle->Phi) < MAX_MATCH_DIST ) {
	    mM++;
	    iMuMatch=iMu;
	  }
	}
	if (iMuMatch!=-1 && iMuMatch<branchMuon->GetEntries()) {
	  mu = (Muon*) branchMuon->At(iMuMatch);
	  hEffPtMu->Fill(genParticle->PT, 1);
	  hEffEtaMu->Fill(genParticle->Eta, 1);
	  hResPtMu->Fill(genParticle->PT, (mu->PT-genParticle->PT)/genParticle->PT);
	  hResEtaMu->Fill(genParticle->Eta, (mu->PT-genParticle->PT)/genParticle->PT);
	}
	else {
	  hEffPtMu->Fill(genParticle->PT, 0);
	  hEffEtaMu->Fill(genParticle->Eta, 0);
	}
      }
      
      if (pid==15) { // TAU
	nT++;
	for (Int_t iJet=0; iJet<branchJet->GetEntries(); iJet++) { //jet loop
	  jet = (Jet*) branchJet->At(iJet);
	  if ( deltaR(jet->Eta, genParticle->Eta, jet->Phi, genParticle->Phi) < MAX_MATCH_DIST ) {
	    mT++;
	    iTauMatch=iJet;
	    if (jet->TauTag>0) mTauTag++;
	  } //end if
	} // end for

	if ( iTauMatch!=-1 && iTauMatch<branchJet->GetEntries() ) { // tau if
	  jet = (Jet*) branchJet->At(iTauMatch);
	  
	  if (jet->TauTag>0) { // tau iff
	    mTauTag++;
	    hEffPtTau->Fill(genParticle->PT, 1);
	    hEffEtaTau->Fill(genParticle->Eta, 1);
	    hResPtTau->Fill(genParticle->PT, (jet->PT-genParticle->PT)/genParticle->PT);
	    hResEtaTau->Fill(genParticle->Eta, (jet->PT-genParticle->PT)/genParticle->PT);
	  } // end tau iff
	  else { // else 
	    hEffPtTau->Fill(genParticle->PT, 0);
	    hEffEtaTau->Fill(genParticle->Eta, 0);
	  } // end else
	  
	} //end tau iff
      } // end TAU if
      
      if (pid==5) {
	nB++;
	for (Int_t iJet=0; iJet<branchJet->GetEntries(); iJet++) { //jet loop
	  jet = (Jet*) branchJet->At(iJet);
	  if ( deltaR(jet->Eta, genParticle->Eta, jet->Phi, genParticle->Phi) < MAX_MATCH_DIST ) {
	    mB++;
	    iBMatch=iJet;
	  }
	}
	if (iBMatch!=-1 && iBMatch<branchJet->GetEntries()) {
	  jet = (Jet*) branchJet->At(iBMatch);
	  if (jet->BTag==1 || jet->BTag==3 || jet->BTag==5 || jet->BTag==7) {
	    mBTag++;
	    hEffPtBJet->Fill(genParticle->PT, 1);
	    hEffEtaBJet->Fill(genParticle->Eta, 1);
	    hResPtBJet->Fill(genParticle->PT, (jet->PT-genParticle->PT)/genParticle->PT);
	    hResEtaBJet->Fill(genParticle->Eta, (jet->PT-genParticle->PT)/genParticle->PT);

	  }
	  
	  else {
	    hEffPtBJet->Fill(genParticle->PT, 0);
	    hEffEtaBJet->Fill(genParticle->Eta, 0);
	  }
	}
      }
    }
  } // end event loop

  TCanvas *c1 = new TCanvas ("c1", "c1", 600, 600);

  hEffPtEle->SetTitle("Ele ID eff as fxn of Pt");
  hEffPtEle->Draw();
  c1->SaveAs("hEffPtEle_phase1.png");
  hEffEtaEle->SetTitle("Ele ID eff as fxn of Eta");
  hEffEtaEle->Draw();
  c1->SaveAs("hEffEtaEle_phase1.png");
  hResPtEle->SetTitle("Ele pt res as fxn of Pt");
  hResPtEle->Draw();
  c1->SaveAs("hResPtEle_phase1.png");
  hResEtaEle->SetTitle("Ele pt res as fxn of Eta");
  hResEtaEle->Draw();
  c1->SaveAs("hResEtaEle_phase1.png");

  hEffPtMu->SetTitle("Mu ID eff as fxn of Pt");
  hEffPtMu->Draw();
  c1->SaveAs("hEffPtMu_phase1.png");
  hEffEtaMu->SetTitle("Mu ID eff as fxn of Eta");
  hEffEtaMu->Draw();
  c1->SaveAs("hEffEtaMu_phase1.png");
  hResPtMu->SetTitle("Mu pt res as fxn of Pt");
  hResPtMu->Draw();
  c1->SaveAs("hResPtMu_phase1.png");
  hResEtaMu->SetTitle("Mu pt res as fxn of Eta");
  hResEtaMu->Draw();
  c1->SaveAs("hResEtaMu_phase1.png");

  hEffPtTau->SetTitle("Tau ID eff as fxn of Pt");
  hEffPtTau->Draw();
  c1->SaveAs("hEffPtTau_phase1.png");
  hEffEtaTau->SetTitle("Tau ID eff as fxn of Eta");
  hEffEtaTau->Draw();
  c1->SaveAs("hEffEtaTau_phase1.png");
  hResPtTau->SetTitle("Tau pt res as fxn of Pt");
  hResPtTau->Draw();
  c1->SaveAs("hResPtTau_phase1.png");
  hResEtaTau->SetTitle("Tau pt res as fxn of Eta");
  hResEtaTau->Draw();
  c1->SaveAs("hResEtaTau_phase1.png");

  hEffPtBJet->SetTitle("BJet ID eff as fxn of Pt");
  hEffPtBJet->Draw();
  c1->SaveAs("hEffPtBJet_phase1.png");
  hEffEtaBJet->SetTitle("BJet ID eff as fxn of Eta");
  hEffEtaBJet->Draw();
  c1->SaveAs("hEffEtaBJet_phase1.png");
  hResPtBJet->SetTitle("BJet pt res as fxn of Pt");
  hResPtBJet->Draw();
  c1->SaveAs("hResPtBJet_phase1.png");
  hResEtaBJet->SetTitle("BJet pt res as fxn of Eta");
  hResEtaBJet->Draw();
  c1->SaveAs("hResEtaBJet_phase1.png");

}

Float_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 ) {
  
  const Float_t pi = 3.14159265358979;
  
  Float_t etaDiff = (eta1-eta2);
  Float_t phiDiff = fabs(phi1-phi2);
  while ( phiDiff>pi ) phiDiff = fabs(phiDiff-2.0*pi);
  
  Float_t deltaRSquared = etaDiff*etaDiff + phiDiff*phiDiff;

  return TMath::Sqrt(deltaRSquared);

}
