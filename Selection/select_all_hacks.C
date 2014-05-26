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

void select_all_hacks(const TString inputfile) {

  TCanvas *c1 = MakeCanvas("c1", "c1", 600, 600);

  const Int_t TAU_ID_CODE = 15;

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
  GenParticle *genEle1=0, *genEle2=0;
  Int_t iGenTau1=-1,    iGenTau2=-1;
  Int_t iGenEle1=-1,    iGenEle2=-1;

  // set up output variables and file
  Float_t allJets=0;
  Float_t nearTaus=0;
  Float_t notNearTaus=0;
  Float_t tauFakes=0;
  Float_t eleFakes=0;
  Float_t eleTot=0;

  Float_t tauTagged=0;
  Float_t bTagged=0;
  Float_t twoTagged=0;

  Float_t xbins[18] = {0,20,30,40,50,60,70,80,90,100,120,140,160,180,200,250,300,400};
  TProfile *tauMis = new TProfile("tauMis", "tauMis", 17, xbins);
  tauMis->SetMarkerColor(kRed);
  tauMis->SetLineColor(kRed);

  for (Int_t iEntry=0; iEntry<numberOfEntries; iEntry++) { // entry loop
    treeReader->ReadEntry(iEntry);

    // ********************
    // RESET
    // ********************
    iGenTau1=-1;    iGenTau2=-1;    
    iGenEle1=-1;    iGenEle2=-1;    
    genTau1=0; genTau2=0;
    genEle1=0; genEle2=0;

    // ********************
    // GEN TAUS
    // ********************
    for (Int_t iParticle=0; iParticle<branchParticle->GetEntries(); iParticle++) { // generator particle loop
      genParticle = (GenParticle*) branchParticle->At(iParticle);
      
      if ( fabs(genParticle->PID) == TAU_ID_CODE ) { // tau switch
	if (!genTau1) {
	  iGenTau1 = iParticle;
	  genTau1 = (GenParticle*) branchParticle->At(iGenTau1);
	}
	else if (!genTau2) {
	  iGenTau2 = iParticle;
	  genTau2 = (GenParticle*) branchParticle->At(iGenTau2);
	}
	else if ((genTau1)&&(genTau2)) {
	  cout << "third tau!" << endl;
	}
      }
      if ( fabs(genParticle->PID) == 11 ) { // ele switch
	if (!genEle1) {
	  iGenEle1 = iParticle;
	  genEle1 = (GenParticle*) branchParticle->At(iGenEle1);
	}
	else if (!genEle2) {
	  iGenEle2 = iParticle;
	  genEle2 = (GenParticle*) branchParticle->At(iGenEle2);
	}
	else if ((genEle1)&&(genEle2)) {
	  cout << "third ele!" << endl;
	}
      }
    } // end gen particle loop
    
    if (genEle1) {
      eleTot+=1;
    }
    if (genEle2) {
      eleTot+=1;
    }

    // ********************
    // FIND FAKES
    // ********************
      
    for (Int_t iJet=0; iJet<branchJet->GetEntries(); iJet++) { // reconstructed jet loop
      jet = (Jet*) branchJet->At(iJet);
      allJets+=1;

      if ((jet->TauTag!=0)&&(jet->BTag!=0)) twoTagged+=1;
      if (jet->TauTag!=0) tauTagged+=1;
      if (jet->BTag!=0) bTagged+=1;
      
      if ((genTau1)&&(deltaR(jet->Eta, genTau1->Eta, jet->Phi, genTau1->Phi) < MAX_MATCH_DIST)) {
	continue;
      }
      if ((genTau2)&&(deltaR(jet->Eta, genTau2->Eta, jet->Phi, genTau2->Phi) < MAX_MATCH_DIST)) {
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

      if ((genEle1)&&(deltaR(jet->Eta, genEle1->Eta, jet->Phi, genEle1->Phi) < MAX_MATCH_DIST)) {
	eleFakes+=1;
      }
      else if ((genEle2)&&(deltaR(jet->Eta, genEle2->Eta, jet->Phi, genEle2->Phi) < MAX_MATCH_DIST)) {
	eleFakes+=1;
      }
	
    } // end reco jet loop
    
  } // end event loop

  c1->SetLogy();

  tauMis->SetTitle("");
  tauMis->GetXaxis()->SetTitle("p_{T}");
  tauMis->GetYaxis()->SetTitle("Fake Rate");
  tauMis->GetYaxis()->SetRangeUser(1e-5,1);
  tauMis->Draw();

  c1->SaveAs("nom_tau_fake.png");
  //tauEle->Divide(eleTot);
  //tauEle->Draw();

  cout << "Tau tag rate  : " << 100*tauTagged/allJets << "%" << endl;
  cout << "B tag rate    : " << 100*bTagged/allJets << "%" << endl;
  cout << "Tau+b tag rate: " << 100*twoTagged/allJets << "%" << endl;
  cout << "b tag rate for taus: " << 100*twoTagged/tauTagged << "%" << endl;
  cout << "Tau fake rate : " << 100*tauFakes/(notNearTaus) << "% " << endl;
  cout << "Total # ele fakes: " << eleFakes << " of " << eleTot << " electrons " << endl;
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
