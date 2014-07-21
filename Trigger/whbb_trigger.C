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

void whbb_trigger(const TString input="whbb.txt",
		  const TString outputfile="/afs/cern.ch/work/j/jlawhorn/whbb.root") {

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
  TClonesArray *branchGenMET = treeReader->UseBranch("GenMissingET");
  TClonesArray *branchHT = treeReader->UseBranch("ScalarHT");
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");

  GenParticle *part=0, *mama=0; Jet *jet=0;
  MissingET *metObj=0;
  ScalarHT *htObj=0;
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
  UInt_t triggerBits126;
  UInt_t triggerBits180;
  UInt_t triggerBits250;
  UInt_t triggerBits350;
  
  TFile *outfile = new TFile(outputfile, "RECREATE");
  TTree *outtree = new TTree("Events", "Events");
  Trigger::Event data;
  outtree->Branch("Events", &data.isReco, "isReco/i:e1pt/F:e1eta:e1phi:e2pt:e2eta:e2phi:m1pt:m1eta:m1phi:m2pt:m2eta:m2phi:g1pt:g1eta:g1phi:g2pt:g2eta:g2phi:t1pt:t1eta:t1phi:t2pt:t2eta:t2phi:j1pt:j1eta:j1phi:j2pt:j2eta:j2phi:j3pt:j3eta:j3phi:j4pt:j4eta:j4phi:mht:ht:triggerBits126/i:triggerBits180:triggerBits250:triggerBits350");

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

    if (branchGenMET) {
      metObj=(MissingET*) branchGenMET->At(0);
      mht=metObj->MET;
    }
    else mht=3;

    if (branchHT) {
      htObj=(ScalarHT*) branchHT->At(0);
      ht=htObj->HT;
    }
    else ht=5;

    //cout << "----" << endl;
    //cout << "      fUID  i  PID  Status  M1  M2  D1  D2" << endl;
    for (Int_t i=0; i<branchParticle->GetEntries(); i++) {
      part = (GenParticle*) branchParticle->At(i);

      if (part->Status==1) continue;
      
      if (fabs(part->PID) == ELE_ID_CODE) {
	//cout << "electron " << part->Status << endl;
	if (part->PT>e1pt) { 
	  if (deltaR(part->Eta, e1eta, part->Phi, e1phi)>MAX_MATCH_DIST) {
	    e2pt=e1pt;
	    e2eta=e1eta;
	    e2phi=e1phi;
	  }
	  e1pt=part->PT;
	  e1eta=part->Eta;
	  e1phi=part->Phi;
	}
	else if (part->PT>e2pt && deltaR(part->Eta, e1eta, part->Phi, e1phi)>MAX_MATCH_DIST) {
	  e2pt=e1pt;
	  e2eta=e1eta;
	  e2phi=e1phi;
	}
      }
      
      if (fabs(part->PID) == MU_ID_CODE) {
	//cout << "muon " << part->Status << endl;
	if (part->PT>m1pt) { 
	  if (deltaR(part->Eta, m1eta, part->Phi, m1phi)>MAX_MATCH_DIST) {
	    m2pt=m1pt;
	    m2eta=m1eta;
	    m2phi=m1phi;
	  }
	  m1pt=part->PT;
	  m1eta=part->Eta;
	  m1phi=part->Phi;
	}
	else if (part->PT>m2pt && deltaR(part->Eta, m1eta, part->Phi, m1phi)>MAX_MATCH_DIST) {
	  m2pt=m1pt;
	  m2eta=m1eta;
	  m2phi=m1phi;
	}
      }
      
      if (fabs(part->PID) == TAU_ID_CODE) {
	//cout << "tau " << part->Status << endl;
	if (part->PT>t1pt) { 
	  if (deltaR(part->Eta, t1eta, part->Phi, t1phi)>MAX_MATCH_DIST) {
	    t2pt=t1pt;
	    t2eta=t1eta;
	    t2phi=t1phi;
	  }
	  t1pt=part->PT;
	  t1eta=part->Eta;
	  t1phi=part->Phi;
	}
	else if (part->PT>t2pt && deltaR(part->Eta, t1eta, part->Phi, t1phi)>MAX_MATCH_DIST) {
	  t2pt=t1pt;
	  t2eta=t1eta;
	  t2phi=t1phi;
	}
      }
    }
    /*    
    for (Int_t i=0; i<branchParticle->GetEntries(); i++) {
      part = (GenParticle*) branchParticle->At(i);
      
      if (fabs(part->PID)==ELE_ID_CODE || fabs(part->PID)==MU_ID_CODE) {
	
	////cout << "lepton cleaning" << endl;
	
	if (deltaR(part->Eta, t1eta, part->Phi, t1phi)<MAX_MATCH_DIST) {
	  //cout << "match to tau 1" << endl;
	  t1pt=-99; t1eta=-99; t1phi=-99;
	}
	if (deltaR(part->Eta, t2eta, part->Phi, t2phi)<MAX_MATCH_DIST) {
	  //cout << "match to tau 2" << endl;
	  t2pt=-99; t2eta=-99; t2phi=-99;
	}
      }
      }*/
    /*
    cout << m1pt << " " << m1eta << " " << m1phi << endl;
    cout << m2pt << " " << m2eta << " " << m2phi << endl;
    cout << e1pt << " " << e1eta << " " << e1phi << endl;
    cout << e2pt << " " << e2eta << " " << e2phi << endl;
    cout << t1pt << " " << t1eta << " " << t1phi << endl;
    cout << t2pt << " " << t2eta << " " << t2phi << endl;
    */
    if (deltaR(e1eta, t1eta, e1phi, t1phi)<MAX_MATCH_DIST) {
      //cout << "e1 matches t1, removing t1" << endl;
      t1pt=-99; t1eta=-99; t1phi=-99;
    }
    if (deltaR(e2eta, t1eta, e2phi, t1phi)<MAX_MATCH_DIST) {
      //cout << "e2 matches t1, removing t1" << endl;
      t1pt=-99; t1eta=-99; t1phi=-99;
    }
    if (deltaR(e1eta, t2eta, e1phi, t2phi)<MAX_MATCH_DIST) {
      //cout << "e1 matches t2, removing t2" << endl;
      t2pt=-99; t2eta=-99; t2phi=-99;
    }
    if (deltaR(e2eta, t2eta, e2phi, t2phi)<MAX_MATCH_DIST) {
      //cout << "e2 matches t2, removing t2" << endl;
      t2pt=-99; t2eta=-99; t2phi=-99;
    }
    if (deltaR(m1eta, t1eta, m1phi, t1phi)<MAX_MATCH_DIST) {
      //cout << "m1 matches t1, removing t1" << endl;
      t1pt=-99; t1eta=-99; t1phi=-99;
    }
    if (deltaR(m2eta, t1eta, m2phi, t1phi)<MAX_MATCH_DIST) {
      //cout << "m2 matches t1, removing t1" << endl;
      t1pt=-99; t1eta=-99; t1phi=-99;
    }
    if (deltaR(m1eta, t2eta, m1phi, t2phi)<MAX_MATCH_DIST) {
      //cout << "m1 matches t2, removing t2" << endl;
      t2pt=-99; t2eta=-99; t2phi=-99;
    }
    if (deltaR(m2eta, t2eta, m2phi, t2phi)<MAX_MATCH_DIST) {
      //cout << "m2 matches t2, removing t2" << endl;
      t2pt=-99; t2eta=-99; t2phi=-99;
    }

    if (t1pt>-99 || t2pt>-99) {
      Float_t t1pt_p=t1pt;
      Float_t t2pt_p=t2pt;

      for (Int_t i=0; i<branchGenJet->GetEntries(); i++) {
	jet = (Jet*) branchGenJet->At(i);

	if (deltaR(jet->Eta, t1eta, jet->Phi, t1phi)<MAX_MATCH_DIST && (t1pt==t1pt_p || jet->PT>t1pt)) {
	  ////cout << "found jet for tau 1" << endl;
	  t1pt=jet->PT;
	  t1eta=jet->Eta;
	  t1phi=jet->Phi;
	}

	if (deltaR(jet->Eta, t2eta, jet->Phi, t2phi)<MAX_MATCH_DIST && (t2pt==t2pt_p || jet->PT>t2pt)) {
	  ////cout << "found jet for tau 2" << endl;
	  t2pt=jet->PT;
	  t2eta=jet->Eta;
	  t2phi=jet->Phi;
	}
	
      }

      if (t1pt_p==t1pt) {
	////cout << "no jet match for tau 1" << endl;
	t1pt=-99; t1eta=-99; t1phi=-99;
      }
      if (t2pt_p==t2pt) {
	////cout << "no jet match for tau 2" << endl;
	t2pt=-99; t2eta=-99; t2phi=-99;
      }

    }
    /*
    cout << m1pt << " " << m1eta << " " << m1phi << endl;
    cout << m2pt << " " << m2eta << " " << m2phi << endl;
    cout << e1pt << " " << e1eta << " " << e1phi << endl;
    cout << e2pt << " " << e2eta << " " << e2phi << endl;
    cout << t1pt << " " << t1eta << " " << t1phi << endl;
    cout << t2pt << " " << t2eta << " " << t2phi << endl;
    */
    for (Int_t i=0; i<branchGenJet->GetEntries(); i++) {
      jet = (Jet*) branchGenJet->At(i);

      if (jet->PT > j1pt) {
	j4pt=j3pt; j4eta=j3eta; j4phi=j3phi;
	j3pt=j2pt; j3eta=j2eta; j3phi=j2phi;
	j2pt=j1pt; j2eta=j1eta; j2phi=j1phi;
	j1pt=jet->PT; j1eta=jet->Eta; j1phi=jet->Phi;
      }
      else if (jet->PT > j2pt) {
	j4pt=j3pt; j4eta=j3eta; j4phi=j3phi;
	j3pt=j2pt; j3eta=j2eta; j3phi=j2phi;
	j2pt=jet->PT; j2eta=jet->Eta; j2phi=jet->Phi;
      }
      else if (jet->PT > j3pt) {
	j4pt=j3pt; j4eta=j3eta; j4phi=j3phi;
	j3pt=jet->PT; j3eta=jet->Eta; j3phi=jet->Phi;
      }
      else if  (jet->PT > j4pt) {
	j4pt=jet->PT; j4eta=jet->Eta; j4phi=jet->Phi;
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
    data.mht=mht;
    data.ht=ht;
    
    triggerBits126 = Trigger::checkTriggers126(data);
    triggerBits180 = Trigger::checkTriggers180(data);
    triggerBits250 = Trigger::checkTriggers250(data);
    triggerBits350 = Trigger::checkTriggers350(data);
    data.triggerBits126=triggerBits126;
    data.triggerBits180=triggerBits180;
    data.triggerBits250=triggerBits250;
    data.triggerBits350=triggerBits350;
    ////cout << "     " << triggerBits << endl;
    outtree->Fill();
    
    //if ( (data.triggerBits & Trigger::kTauMu) >0 ) {
    ////cout << "mutau" << endl;
    //}

    TString tt="t1pt>0 && t2pt>0";
    TString mt="m1pt>0 && (t1pt>0 || t2pt>0)";
    TString et="e1pt>0 && (t1pt>0 || t2pt>0)";
    TString em="m1pt>0 && e1pt>0";
    TString ee="e1pt>0 && e2pt>0";
    TString mm="m1pt>0 && m2pt>0";

    /*    if ( m1pt>0 && (t1pt>0||t2pt>0) && e1pt>0) {
      cout << "problem child" << endl;
    }
    if (!(t1pt>0 && t2pt>0)*!(m1pt>0 && (t1pt>0 || t2pt>0))*!(e1pt>0 && (t1pt>0 || t2pt>0))*!(m1pt>0 && e1pt>0)*!(e1pt>0 && e2pt>0)*!(m1pt>0 && m2pt>0) ) {
      cout << "!!!!" << endl;
      cout << "      fUID  i  PID  Status  M1  M2  D1  D2" << endl;
      for (Int_t i=0; i<branchParticle->GetEntries(); i++) {
	part = (GenParticle*) branchParticle->At(i);
	//if ( part->M1>20 ) continue;

	cout << part->fUniqueID << "  " << setw(3) << i << "  " << setw(3) << part->PID << "  " << setw(3) << part->Status << "  " << setw(3) << part->M1 << "  " << setw(3) << part->M2 << "  " << setw(3) << part->D1 << "  " << setw(3) << part->D2 << endl;
	}
      
      }*/
    
    
  }
  
  outfile->Write();
  outfile->Close();
  delete outfile;

  
}
