// *****
// For efficiency studies
// *****
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TVector2.h>
#include <TMath.h>
#include <TChain.h>
#include <TProfile.h>
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
#include "/afs/cern.ch/user/k/klawhorn/DelphesDiHiggs/Utils/HttStyles.h"
#include "/afs/cern.ch/user/k/klawhorn/DelphesDiHiggs/Utils/CPlot.hh"
#endif

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;

Float_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 );

void select_all(const TString inputfile, const Int_t count=0) {

  TCanvas *c1 = MakeCanvas("c1", "c1", 600, 600);

  //cout << count << endl;

  const Int_t TAU_ID_CODE = 15;
  const Int_t B_ID_CODE = 5;

  const Float_t MAX_MATCH_DIST = 0.5;

  // event types
  enum { HH=0, TT, ZH, WH, WW, ZZ, ZW, ETC };

  // tau decay modes
  enum { hadron=1, electron, muon };

  // read input input file
  TChain chain("Delphes");

  ifstream ifs;
  ifs.open(inputfile.Data()); assert(ifs.is_open());
  string line;
  while(getline(ifs,line)) {
    string fname;
    stringstream ss(line);
    ss >> fname;
    chain.Add(TString(fname));
  }
  ifs.close();

  //chain.Add(inputfile);

  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  //TClonesArray *branchEvent = treeReader->UseBranch("Event");

  //set up loop variables
  GenParticle *genParticle=0;
  Jet *jet=0;
  // set up storage variables
  GenParticle *genTau1=0, *genTau2=0;
  GenParticle *genB1=0, *genB2=0;
  Int_t iGenTau1=-1,    iGenTau2=-1;
  Int_t iGenB1=-1,    iGenB2=-1;
  Jet *jetTau1=0, *jetTau2=0;
  Int_t iJetTau1=-1, iJetTau2=-1;
  Jet *jetB1=0, *jetB2=0;
  Int_t iJetB1=-1, iJetB2=-1;

  // set up output variables and file
  Float_t allJets=0;
  Float_t nearTaus=0;
  Float_t nearTausAndTauTagged=0;
  Float_t nearTausAndNotTauTagged=0;
  Float_t notNearTaus=0;

  Float_t nearBs=0;
  Float_t nearBsAndBTagged=0;
  Float_t nearBsAndNotBTagged=0;
  Float_t notNearBs=0;

  Float_t tauFakes=0;
  Float_t bFakes=0;
  Float_t eleFakes=0;
  Float_t eleTot=0;
  Float_t totalTaus=0;
  Float_t totalBs=0;

  Float_t tauTagged=0;
  Float_t bTagged=0;
  Float_t twoTagged=0;

  char outfile[50];
  sprintf(outfile, "/afs/cern.ch/work/k/klawhorn/test/hist_%i.root", count);

  //cout << outfile << endl;

  //TFile *file = new TFile(outfile, "RECREATE");

  Float_t xbins[18] = {0,20,30,40,50,60,70,80,90,100,120,140,160,180,200,250,300,400};

  TProfile *tauMis = new TProfile("tauMis", "tauMis", 17, xbins);
  tauMis->SetMarkerColor(kRed);
  tauMis->SetLineColor(kRed);

  TProfile *tauJet = new TProfile("tauJet", "tauJet", 17, xbins);
  tauJet->SetMarkerColor(kRed);
  tauJet->SetLineColor(kRed);

  TProfile *tauRes = new TProfile("tauRes", "tauRes", 17, xbins);
  tauRes->SetMarkerColor(kRed);
  tauRes->SetLineColor(kRed);

  TProfile *bMis = new TProfile("bMis", "bMis", 17, xbins);
  bMis->SetMarkerColor(kRed);
  bMis->SetLineColor(kRed);

  TProfile *bJet = new TProfile("bJet", "bJet", 17, xbins);
  bJet->SetMarkerColor(kRed);
  bJet->SetLineColor(kRed);

  TProfile *bRes = new TProfile("bRes", "bRes", 17, xbins);
  bRes->SetMarkerColor(kRed);
  bRes->SetLineColor(kRed);

  for (Int_t iEntry=0; iEntry<numberOfEntries; iEntry++) { // entry loop
    treeReader->ReadEntry(iEntry);

    // ********************
    // RESET
    // ********************
    iGenTau1=-1;    iGenTau2=-1;    
    genTau1=0; genTau2=0;
    iJetTau1=-1;    iJetTau2=-1;    
    jetTau1=0; jetTau2=0;

    iGenB1=-1;    iGenB2=-1;    
    genB1=0; genB2=0;
    iJetB1=-1;    iJetB2=-1;    
    jetB1=0; jetB2=0;

    // ********************
    // GEN TAUS
    // ********************
    for (Int_t iParticle=0; iParticle<branchParticle->GetEntries(); iParticle++) { // generator particle loop
      genParticle = (GenParticle*) branchParticle->At(iParticle);
      
      if ( fabs(genParticle->PID) == TAU_ID_CODE ) { // tau switch
	if (!genTau1) {
	  iGenTau1 = iParticle;
	  genTau1 = (GenParticle*) branchParticle->At(iGenTau1);
	  totalTaus+=1;
	}
	else if (!genTau2) {
	  iGenTau2 = iParticle;
	  genTau2 = (GenParticle*) branchParticle->At(iGenTau2);
	  totalTaus+=1;
	}
	else if ((genTau1)&&(genTau2)) {
	  cout << "third tau!" << endl;
	}
      }
      if ( fabs(genParticle->PID) == B_ID_CODE ) { // b switch
	if (!genB1) {
	  iGenB1 = iParticle;
	  genB1 = (GenParticle*) branchParticle->At(iGenB1);
	  totalBs+=1;
	}
	else if (!genB2) {
	  iGenB2 = iParticle;
	  genB2 = (GenParticle*) branchParticle->At(iGenB2);
	  totalBs+=1;
	}
	else if ((genB1)&&(genB2)) {
	  cout << "third b!" << endl;
	}
      }

    } // end gen particle loop
    
    // ********************
    // TAUS
    // ********************
      
    for (Int_t iJet=0; iJet<branchJet->GetEntries(); iJet++) { // reconstructed jet loop
      jet = (Jet*) branchJet->At(iJet);
      allJets+=1;

      if ((jet->TauTag!=0)&&(jet->BTag!=0)) twoTagged+=1;
      if (jet->TauTag!=0) tauTagged+=1;
      if (jet->BTag!=0) bTagged+=1;
      
      if ((genTau1)&&(deltaR(jet->Eta, genTau1->Eta, jet->Phi, genTau1->Phi) < MAX_MATCH_DIST)) {
	iJetTau1=iJet;
	jetTau1 = (Jet*) branchJet->At(iJetTau1);
	nearTaus+=1;
	if (jet->TauTag!=0) {
	  nearTausAndTauTagged+=1;
	  tauJet->Fill(genTau1->PT, 1);
	  tauRes->Fill(genTau1->PT, (jetTau1->PT-genTau1->PT)/genTau1->PT);
	}
	else {
	  nearTausAndNotTauTagged+=1;
	  tauJet->Fill(genTau1->PT, 0);
	}
	continue;
      }
      if ((genTau2)&&(deltaR(jet->Eta, genTau2->Eta, jet->Phi, genTau2->Phi) < MAX_MATCH_DIST)) {
	iJetTau2=iJet;
	jetTau2 = (Jet*) branchJet->At(iJetTau2);
	nearTaus+=1;
	if (jet->TauTag!=0) {
	  nearTausAndTauTagged+=1; 
	  tauJet->Fill(genTau2->PT, 1);
	  tauRes->Fill(genTau2->PT, (jetTau2->PT-genTau2->PT)/genTau2->PT);
	}
	else {
	  nearTausAndNotTauTagged+=1;
	  tauJet->Fill(genTau2->PT, 0);
	}
	continue;
      }
      notNearTaus+=1;
      if (jet->TauTag!=0) {
	tauFakes+=1;
	tauMis->Fill(jet->PT, 1);
      }
      else {
	tauMis->Fill(jet->PT, 0);
      }
    } // end reco jet loop

    // ********************
    // Bs
    // ********************

    for (Int_t iJet=0; iJet<branchJet->GetEntries(); iJet++) { // reconstructed jet loop
      jet = (Jet*) branchJet->At(iJet);

      if ((genB1)&&(deltaR(jet->Eta, genB1->Eta, jet->Phi, genB1->Phi) < MAX_MATCH_DIST)) {
	iJetB1=iJet;
	jetB1 = (Jet*) branchJet->At(iJetB1);
	nearBs+=1;
	if (jet->BTag!=0) {
	  nearBsAndBTagged+=1;
	  bJet->Fill(genB1->PT, 1);
	  bRes->Fill(genB1->PT, (jetB1->PT-genB1->PT)/genB1->PT);
	}
	else {
	  nearBsAndNotBTagged+=1;
	  bJet->Fill(genB1->PT, 0);
	}
	continue;
      }
      if ((genB2)&&(deltaR(jet->Eta, genB2->Eta, jet->Phi, genB2->Phi) < MAX_MATCH_DIST)) {
	iJetB2=iJet;
	jetB2 = (Jet*) branchJet->At(iJetB2);
	nearBs+=1;
	if (jet->BTag!=0) {
	  nearBsAndBTagged+=1; 
	  bJet->Fill(genB2->PT, 1);
	  bRes->Fill(genB2->PT, (jetB2->PT-genB2->PT)/genB2->PT);
	}
	else {
	  nearBsAndNotBTagged+=1;
	  bJet->Fill(genB2->PT, 0);
	}
	continue;
      }
      notNearBs+=1;
      if (jet->BTag!=0) {
	bFakes+=1;
	bMis->Fill(jet->PT, 1);
      }
      else {
	bMis->Fill(jet->PT, 0);
      }
    } // end reco jet loop
    
  } // end event loop

  tauMis->SetTitle("");
  tauMis->GetXaxis()->SetTitle("Reco p_{T}");
  tauMis->GetYaxis()->SetTitle("#tau-tag Fake Rate");
  tauMis->Draw();
  c1->SaveAs("tau_tag_fake.png");

  tauJet->SetTitle("");
  tauJet->GetXaxis()->SetTitle("Generator p_{T}");
  tauJet->GetYaxis()->SetTitle("#tau-tag Efficiency");
  tauJet->Draw();
  c1->SaveAs("tau_tag_eff.png");

  tauRes->SetTitle("");
  tauRes->GetXaxis()->SetTitle("Generator p_{T}");
  tauRes->GetYaxis()->SetTitle("Resolution (reco-gen)/(gen)");
  tauRes->Draw();
  c1->SaveAs("tau_tag_res.png");

  bMis->SetTitle("");
  bMis->GetXaxis()->SetTitle("Reco p_{T}");
  bMis->GetYaxis()->SetTitle("b-tag Fake Rate");
  bMis->Draw();
  c1->SaveAs("b_tag_fake.png");

  bJet->SetTitle("");
  bJet->GetXaxis()->SetTitle("Generator p_{T}");
  bJet->GetYaxis()->SetTitle("b-tag Efficiency");
  bJet->Draw();
  c1->SaveAs("b_tag_eff.png");

  bRes->SetTitle("");
  bRes->GetXaxis()->SetTitle("Generator p_{T}");
  bRes->GetYaxis()->SetTitle("Resolution (reco-gen)/(gen)");
  bRes->Draw();
  c1->SaveAs("b_tag_res.png");

  cout << "Total taus: " << totalTaus << endl;
  cout << "N events  : " << numberOfEntries << endl;
  cout << "Reco jets near taus: " << nearTaus << endl;
  cout << "  and tau tagged: " << nearTausAndTauTagged << endl;
  cout << "  and not tau tagged: " << nearTausAndNotTauTagged << endl;
  cout << "Tau efficiency: " << nearTausAndTauTagged/nearTaus << endl;

  cout << "Total bs: " << totalBs << endl;
  cout << "N events  : " << numberOfEntries << endl;
  cout << "Reco jets near bs: " << nearBs << endl;
  cout << "  and b tagged: " << nearBsAndBTagged << endl;
  cout << "  and not b tagged: " << nearBsAndNotBTagged << endl;
  cout << "B efficiency: " << nearBsAndBTagged/nearBs << endl;

  cout << "Tau tag rate  : " << 100*tauTagged/allJets << "%" << endl;
  cout << "B tag rate    : " << 100*bTagged/allJets << "%" << endl;
  cout << "Tau+b tag rate: " << 100*twoTagged/allJets << "%" << endl;
  cout << endl;

}

Float_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 ) {

  const Float_t pi = 3.14159265358979;

  Float_t etaDiff = (eta1-eta2);
  Float_t phiDiff = fabs(phi1-phi2);
  while ( phiDiff>pi ) phiDiff = fabs(phiDiff-2.0*pi);

  Float_t deltaRSquared = etaDiff*etaDiff + phiDiff*phiDiff;

  return TMath::Sqrt(deltaRSquared);

}
