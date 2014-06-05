#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
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
#include "modules/Delphes.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "classes/DelphesClasses.h"
#endif

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;

Float_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 );

void select_bare(const TString inputfile="qcd_phaseII_conf4.txt", TString outputfile="qcd.root")
{

  Float_t maxR = 0.4;
  //Float_t maxR = 0.25;

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
  TClonesArray *branchGenJet = treeReader->UseBranch("GenJet");
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");

  TClonesArray *branchMissingET = treeReader->UseBranch("MissingET");
  TClonesArray *branchGenMissingET = treeReader->UseBranch("GenMissingET");
  TClonesArray *branchPileUpJetIDMissingET = treeReader->UseBranch("PileUpJetIDMissingET");

  Int_t isPU=1;
  Int_t genPass=1;

  //set up loop variables
  Jet *genJet=0;
  Jet *jet=0;
  GenParticle *part=0;

  MissingET *met1=0;
  MissingET *met2=0;
  MissingET *met3=0;
  
  Float_t pt=0, eta=0, phi=0;
  Float_t beta=0, betastar=0, mdrsq=0;
  Float_t met=0, genmet=0, puidmet=0;
  Float_t metphi=0, genmetphi=0, puidmetphi=0;
  Int_t pu=0;

  TFile *outfile = new TFile(outputfile, "RECREATE");
  TTree *outtree = new TTree("Events", "Events");

  outtree->Branch("pt", &pt, "pt/f");
  outtree->Branch("eta", &eta, "eta/f");
  outtree->Branch("phi", &phi, "phi/f");
  outtree->Branch("beta", &beta, "beta/f");
  outtree->Branch("betastar", &betastar, "betastar/f");
  outtree->Branch("mdrsq", &mdrsq, "mdrsq/f");
  outtree->Branch("met", &met, "met/f");
  outtree->Branch("genmet", &genmet, "genmet/f");
  outtree->Branch("puidmet", &puidmet, "puidmet/f");
  outtree->Branch("metphi", &metphi, "metphi/f");
  outtree->Branch("genmetphi", &genmetphi, "genmetphi/f");
  outtree->Branch("puidmetphi", &puidmetphi, "puidmetphi/f");
  outtree->Branch("pu", &pu, "pu/i");
  
  Int_t skip=0;

  for (Int_t iEntry=0; iEntry<numberOfEntries; iEntry++) { // entry loop
    treeReader->ReadEntry(iEntry);

    pt=0; eta=0; phi=0; beta=0; betastar=0; mdrsq=0; pu=0;
    met=0; genmet=0; puidmet=0; metphi=0; genmetphi=0; puidmetphi=0; 
    
    for (Int_t iJet=0; iJet<branchJet->GetEntries(); iJet++) { // reconstructed jet loop
      jet = (Jet*) branchJet->At(iJet);

      if (jet->PT < 20) continue;

      skip=0;

      isPU=1;
      genPass=0;

      for (Int_t iPart=0; iPart<branchParticle->GetEntries(); iPart++) {
	part = (GenParticle*) branchParticle->At(iPart);

	if (fabs(part->PID) != 15) continue;

	if (deltaR(jet->Eta, part->Eta, jet->Phi, part->Phi)<maxR) {
	  skip=1;
	}
      }

      if (skip==1) continue;

      for (Int_t iGenJet=0; iGenJet<branchGenJet->GetEntries(); iGenJet++) { // generator level jet loop
	genJet = (Jet*) branchGenJet->At(iGenJet);
	if (deltaR(jet->Eta, genJet->Eta, jet->Phi, genJet->Phi)<maxR) { 
	  if (isPU==0) {
	    //cout << "found second gen jet?" << endl;
	    if ((genPass!=1) && (genJet->PT>10)) genPass=1;
	  }
	  else if (genJet->PT>10) {
	    genPass=1;
	    isPU=0; 
	  }
	}
      }

      pt=jet->PT;
      eta=jet->Eta;
      phi=jet->Phi;
      beta=jet->Beta;
      betastar=jet->BetaStar;
      mdrsq=jet->MeanSqDeltaR;
      pu=isPU;

      met1 = (MissingET*) branchMissingET->At(0);
      met2 = (MissingET*) branchGenMissingET->At(0);
      met3 = (MissingET*) branchPileUpJetIDMissingET->At(0);

      met=met1->MET;
      metphi=met1->Phi;
      genmet=met2->MET;
      genmetphi=met2->Phi;
      puidmet=met3->MET;
      puidmetphi=met3->Phi;

      outtree->Fill();

    } // end reco jet loop

  } // end event loop

  outfile->Write();
  outfile->Close();

}

Float_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 ) {

  const Float_t pi = 3.14159265358979;

  Float_t etaDiff = (eta1-eta2);
  Float_t phiDiff = fabs(phi1-phi2);
  while ( phiDiff>pi ) phiDiff = fabs(phiDiff-2.0*pi);

  Float_t deltaRSquared = etaDiff*etaDiff + phiDiff*phiDiff;

  return TMath::Sqrt(deltaRSquared);

}
