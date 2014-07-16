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
#include <sstream>
#include "Math/LorentzVector.h"
#include "TLorentzVector.h"

#include "modules/Delphes.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "classes/DelphesClasses.h"

#include "triggerMenu.h"

#endif

void generalTrigger(const TString input="htt.txt",
		    const TString outputfile="/afs/cern.ch/work/j/jlawhorn/htt.root") {

  // declare constants
  //const Double_t MUON_MASS = 0.105658369;
  //const Double_t ELE_MASS  = 0.000511;
  //const Double_t TAU_MASS  = 1.77682;
  
  const Int_t GAM_ID_CODE = 22;
  const Int_t TAU_ID_CODE = 15;
  const Int_t MU_ID_CODE = 13;
  const Int_t ELE_ID_CODE = 11;
  
  const Float_t MAX_MATCH_DIST = 0.4;
  
  TChain chain("Delphes");
  
  ifstream ifs;
  ifs.open(input.Data()); assert(ifs.is_open());
  string line;
  while(getline(ifs,line)) {
    TString inputfile;
    stringstream ss(line);
    ss >> inputfile;
    chain.Add(inputfile);
  }
  ifs.close();
  
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t n = treeReader->GetEntries();
  
  //TClonesArray *branchJet = treeReader->UseBranch("Jet");
  
  //if (!(branchJet)) {
  //cout << "File broken" << endl;
  //return;
  //}
  
  //TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  //TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchGenJet = treeReader->UseBranch("GenJet");
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");

  GenParticle *part=0; Jet *jet=0;
  //Electron *ele=0; Muon *mu=0;

  UInt_t isReco=0;
  Float_t e1pt, e1eta, e1phi;
  Float_t e2pt, e2eta, e2phi;
  Float_t m1pt, m1eta, m1phi;
  Float_t m2pt, m2eta, m2phi;
  Float_t g1pt, g1eta, g1phi;
  Float_t g2pt, g2eta, g2phi;
  Float_t t1pt, t1eta, t1phi;
  Float_t t2pt, t2eta, t2phi;
  Float_t j1pt, j1eta, j1phi;
  Float_t j2pt, j2eta, j2phi;
  Float_t j3pt, j3eta, j3phi;
  Float_t j4pt, j4eta, j4phi;
  Float_t mht, ht;
  UInt_t triggerBits;
  
  TFile *outfile = new TFile(outputfile, "RECREATE");
  TTree *outtree = new TTree("Events", "Events");
  Trigger::Event data;
  outtree->Branch("Events", &data.isReco, "isReco/i:e1pt/F:e1eta:e1phi:e2pt:e2eta:e2phi:m1pt:m1eta:m1phi:m2pt:m2eta:m2phi:g1pt:g1eta:g1phi:g2pt:g2eta:g2phi:t1pt:t1eta:t1phi:t2pt:t2eta:t2phi:j1pt:j1eta:j1phi:j2pt:j2eta:j2phi:j3pt:j3eta:j3phi:j4pt:j4eta:j4phi:mht:ht:triggerBits/i");

  //n=50;
  
  for (Int_t iEntry=0; iEntry<n; iEntry++) { // entry loop                                                                          
    treeReader->ReadEntry(iEntry);

    e1pt=-99; e1eta=-99; e1phi=-99; 
    e2pt=-99; e2eta=-99; e2phi=-99;
    m1pt=-99; m1eta=-99; m1phi=-99; 
    m2pt=-99; m2eta=-99; m2phi=-99; 
    g1pt=-99; g1eta=-99; g1phi=-99; 
    g2pt=-99; g2eta=-99; g2phi=-99; 
    t1pt=-99; t1eta=-99; t1phi=-99; 
    t2pt=-99; t2eta=-99; t2phi=-99; 
    j1pt=-99; j1eta=-99; j1phi=-99; 
    j2pt=-99; j2eta=-99; j2phi=-99; 
    j3pt=-99; j3eta=-99; j3phi=-99; 
    j4pt=-99; j4eta=-99; j4phi=-99; 
    mht=-99; ht=-99; 

    //cout << "----" << endl;

    for (Int_t i=0; i<branchParticle->GetEntries(); i++) {
      part = (GenParticle*) branchParticle->At(i);
      //cout << part->PID << " " << part->Status << " " << part->PT << " " << part->Eta << endl;
      //if (part->Status!=3) continue;

      if (fabs(part->PID) == ELE_ID_CODE && part->PT > e1pt ) {
	if (e1pt>e2pt) {
	  e2pt=e1pt;
	  e2eta=e1eta;
	  e2phi=e1phi;
	}
	e1pt=part->PT;
	e1eta=part->Eta;
	e1phi=part->Phi;
      }
      else if (fabs(part->PID) == ELE_ID_CODE && part->PT > e2pt) {
	e2pt=part->PT;
	e2eta=part->Eta;
	e2eta=part->Phi;
      }
      else if (fabs(part->PID) == MU_ID_CODE && part->PT > m1pt ) {
	if (m1pt>m2pt) {
	  m2pt=m1pt;
	  m2eta=m1eta;
	  m2phi=m1phi;
	}
	m1pt=part->PT;
	m1eta=part->Eta;
	m1phi=part->Phi;
      }
      else if (fabs(part->PID) == MU_ID_CODE && part->PT > m2pt) {
	m2pt=part->PT;
	m2eta=part->Eta;
	m2phi=part->Phi;
      }
      else if (fabs(part->PID) == TAU_ID_CODE && part->PT > t1pt ) {
	if (t1pt>t2pt) {
	  t2pt=t1pt;
	  t2eta=t1eta;
	  t2phi=t1phi;
	}
	t1pt=part->PT;
	t1eta=part->Eta;
	t1phi=part->Phi;
      }
      else if (fabs(part->PID) == TAU_ID_CODE && part->PT > t2pt) {
	t2pt=part->PT;
	t2eta=part->Eta;
	t2phi=part->Phi;
      }
      else if (fabs(part->PID) == GAM_ID_CODE && part->PT > g1pt ) {
	if (g1pt>g2pt) {
	  g2pt=g1pt;
	  g2eta=g1eta;
	  g2phi=g1phi;
	}
	g1pt=part->PT;
	g1eta=part->Eta;
	g1phi=part->Phi;
      }
      else if (fabs(part->PID) == GAM_ID_CODE && part->PT > g2pt) {
	g2pt=part->PT;
	g2eta=part->Eta;
	g2phi=part->Phi;
      }
    }
    
    for (Int_t i=0; i<branchParticle->GetEntries(); i++) {
      
      if (fabs(part->PID)!=ELE_ID_CODE && fabs(part->PID)!=MU_ID_CODE) continue;
      
      if (deltaR(part->Eta, t1eta, part->Phi, t1phi)<MAX_MATCH_DIST) {
	t1pt=-99; t1eta=-99; t1phi=-99;
      }
      if (deltaR(part->Eta, t2eta, part->Phi, t2phi)<MAX_MATCH_DIST) {
	t2pt=-99; t2eta=-99; t2phi=-99;
      }
    }

    if (t1pt>-99 || t2pt>-99) {
      Float_t t1pt_p=t1pt;
      Float_t t2pt_p=t2pt;

      for (Int_t i=0; i<branchGenJet->GetEntries(); i++) {
	jet = (Jet*) branchGenJet->At(i);

	if (deltaR(jet->Eta, t1eta, jet->Phi, t1phi)<MAX_MATCH_DIST && (t1pt==t1pt_p || jet->PT>t1pt)) {
	  t1pt=jet->PT;
	  t1eta=jet->Eta;
	  t1phi=jet->Phi;
	}

	if (deltaR(jet->Eta, t2eta, jet->Phi, t2phi)<MAX_MATCH_DIST && (t2pt==t2pt_p || jet->PT>t2pt)) {
	  t2pt=jet->PT;
	  t2eta=jet->Eta;
	  t2phi=jet->Phi;
	}
	
      }
      
    }
    
    data.isReco=0;
    data.e1pt=e1pt;
    data.e1eta=e1eta;
    data.e1phi=e1phi;
    data.e2pt=e2pt;
    data.e2eta=e2eta;
    data.e2phi=e2phi;
    data.m1pt=m1pt;
    data.m1eta=m1eta;
    data.m1phi=m1phi;
    data.m2pt=m2pt;
    data.m2eta=m2eta;
    data.m2phi=m2phi;
    data.g1pt=g1pt;
    data.g1eta=g1eta;
    data.g1phi=g1phi;
    data.g2pt=g2pt;
    data.g2eta=g2eta;
    data.g2phi=g2phi;
    data.t1pt=t1pt;
    data.t1eta=t1eta;
    data.t1phi=t1phi;
    data.t2pt=t2pt;
    data.t2eta=t2eta;
    data.t2phi=t2phi;
    data.j1pt=j1pt;
    data.j1eta=j1eta;
    data.j1phi=j1phi;
    data.j2pt=j2pt;
    data.j2eta=j2eta;
    data.j2phi=j2phi;
    data.j3pt=j3pt;
    data.j3eta=j3eta;
    data.j3phi=j3phi;
    data.j4pt=j4pt;
    data.j4eta=j4eta;
    data.j4phi=j4phi;
    data.mht=7;
    data.ht=234;
    
    triggerBits = Trigger::checkTriggers(data);
    data.triggerBits=triggerBits;
    //cout << "     " << triggerBits << endl;
    outtree->Fill();

    if ( (data.triggerBits & Trigger::kTauMu) >0 ) {
      cout << "mutau" << endl;
    }
  }
  
  outfile->Write();
  outfile->Close();
  delete outfile;

  
}
