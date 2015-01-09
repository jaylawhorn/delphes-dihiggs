#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TVector2.h>
#include <TMath.h>
#include <TChain.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "../Utils/hhMVA.h"

#endif

Float_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 );

void flag_events(const TString input="/afs/cern.ch/work/j/jlawhorn/HHToTTBB_14TeV.root",
		 //const TString input="/afs/cern.ch/work/j/jlawhorn/tt.root",
		 const TString output="/afs/cern.ch/work/j/jlawhorn/HHToTTBB-small.root") {
		 //const TString output="/afs/cern.ch/work/j/jlawhorn/tt-small.root") {

  // read input file
  TFile *file = new TFile(input, "READ");
  TTree *tree = (TTree *) file->Get("Events");

  UInt_t isBBTT;

  UInt_t tauCat1=0, tauCat2=0;
  UInt_t bTag1=0, bTag2=0;

  Float_t ptTau1, ptTau2, ptB1, ptB2;
  Float_t etaTau1, etaTau2, etaB1, etaB2;
  Float_t phiTau1, phiTau2, phiB1, phiB2;
  Float_t mTau1, mTau2, mB1, mB2;

  Float_t ptTau1_gen, ptTau2_gen, ptB1_gen, ptB2_gen, ptB3_gen, ptB4_gen;
  Float_t etaTau1_gen, etaTau2_gen, etaB1_gen, etaB2_gen, etaB3_gen, etaB4_gen;
  Float_t phiTau1_gen, phiTau2_gen, phiB1_gen, phiB2_gen, phiB3_gen, phiB4_gen;
  Float_t mTau1_gen, mTau2_gen, mB1_gen, mB2_gen, mB3_gen, mB4_gen;

  Float_t ptTau1_genJet, ptTau2_genJet, etaTau1_genJet, etaTau2_genJet, phiTau1_genJet, phiTau2_genJet, mTau1_genJet, mTau2_genJet;

  TFile *outFile = new TFile(output, "RECREATE");

  tree->SetBranchAddress("isBBTT",         &isBBTT);

  tree->SetBranchAddress("ptTau1",         &ptTau1);
  tree->SetBranchAddress("etaTau1",        &etaTau1);
  tree->SetBranchAddress("phiTau1",        &phiTau1);
  tree->SetBranchAddress("mTau1",          &mTau1);
  tree->SetBranchAddress("tauCat1",        &tauCat1);

  tree->SetBranchAddress("ptTau2",         &ptTau2);
  tree->SetBranchAddress("etaTau2",        &etaTau2);
  tree->SetBranchAddress("phiTau2",        &phiTau2);
  tree->SetBranchAddress("mTau2",          &mTau2);
  tree->SetBranchAddress("tauCat2",        &tauCat2);

  tree->SetBranchAddress("ptB1",           &ptB1);
  tree->SetBranchAddress("etaB1",          &etaB1);
  tree->SetBranchAddress("phiB1",          &phiB1);
  tree->SetBranchAddress("mB1",            &mB1);
  tree->SetBranchAddress("bTag1",          &bTag1);

  tree->SetBranchAddress("ptB2",           &ptB2);
  tree->SetBranchAddress("etaB2",          &etaB2);
  tree->SetBranchAddress("phiB2",          &phiB2);
  tree->SetBranchAddress("mB2",            &mB2);
  tree->SetBranchAddress("bTag2",          &bTag2);

  tree->SetBranchAddress("ptTau1_gen",     &ptTau1_gen);
  tree->SetBranchAddress("etaTau1_gen",    &etaTau1_gen);
  tree->SetBranchAddress("phiTau1_gen",    &phiTau1_gen);
  tree->SetBranchAddress("mTau1_gen",      &mTau1_gen);

  tree->SetBranchAddress("ptTau2_gen",     &ptTau2_gen);
  tree->SetBranchAddress("etaTau2_gen",    &etaTau2_gen);
  tree->SetBranchAddress("phiTau2_gen",    &phiTau2_gen);
  tree->SetBranchAddress("mTau2_gen",      &mTau2_gen);

  tree->SetBranchAddress("ptTau1_genJet",  &ptTau1_genJet);
  tree->SetBranchAddress("etaTau1_genJet", &etaTau1_genJet);
  tree->SetBranchAddress("phiTau1_genJet", &phiTau1_genJet);
  tree->SetBranchAddress("mTau1_genJet",   &mTau1_genJet);

  tree->SetBranchAddress("ptTau2_genJet",  &ptTau2_genJet);
  tree->SetBranchAddress("etaTau2_genJet", &etaTau2_genJet);
  tree->SetBranchAddress("phiTau2_genJet", &phiTau2_genJet);
  tree->SetBranchAddress("mTau2_genJet",   &mTau2_genJet);

  tree->SetBranchAddress("ptB1_gen",       &ptB1_gen);
  tree->SetBranchAddress("etaB1_gen",      &etaB1_gen);
  tree->SetBranchAddress("phiB1_gen",      &phiB1_gen);
  tree->SetBranchAddress("mB1_gen",        &mB1_gen);

  tree->SetBranchAddress("ptB2_gen",       &ptB2_gen);
  tree->SetBranchAddress("etaB2_gen",      &etaB2_gen);
  tree->SetBranchAddress("phiB2_gen",      &phiB2_gen);
  tree->SetBranchAddress("mB2_gen",        &mB2_gen);

  tree->SetBranchAddress("ptB3_gen",       &ptB3_gen);
  tree->SetBranchAddress("etaB3_gen",      &etaB3_gen);
  tree->SetBranchAddress("phiB3_gen",      &phiB3_gen);
  tree->SetBranchAddress("mB3_gen",        &mB3_gen);

  tree->SetBranchAddress("ptB4_gen",       &ptB4_gen);
  tree->SetBranchAddress("etaB4_gen",      &etaB4_gen);
  tree->SetBranchAddress("phiB4_gen",      &phiB4_gen);
  tree->SetBranchAddress("mB4_gen",        &mB4_gen);

  tree->GetEntry(0);

  TTree *outtree=(TTree*)tree->GetTree()->CloneTree(0);

  Int_t rTau1=0, rTau2=0, rB1=0, rB2=0;

  outtree->Branch("rTau1", &rTau1, "rTau1/i");
  outtree->Branch("rTau2", &rTau2, "rTau2/i");
  outtree->Branch("rB1", &rB1, "rB1/i");
  outtree->Branch("rB2", &rB2, "rB2/i");

  for (Int_t i=0; i<tree->GetEntries(); i++) {
    tree->GetEntry(i);

    rTau1=0; rTau2=0;
    rB1=0; rB2=0;
    
    if (isBBTT!=1) continue;

    if (deltaR(etaTau1, etaTau1_gen, phiTau1, phiTau1_gen)<0.4) rTau1=1;
    else if (deltaR(etaTau1, etaTau2_gen, phiTau1, phiTau2_gen)<0.4) rTau1=1;

    if (deltaR(etaTau2, etaTau2_gen, phiTau2, phiTau2_gen)<0.4) rTau2=1;
    else if (deltaR(etaTau2, etaTau1_gen, phiTau2, phiTau1_gen)<0.4) rTau2=1;
    
    if (deltaR(etaB1, etaB1_gen, phiB1, phiB1_gen)<0.4) rB1=1;
    else if (deltaR(etaB1, etaB2_gen, phiB1, phiB2_gen)<0.4) rB1=1;
    else if (deltaR(etaB1, etaB3_gen, phiB1, phiB3_gen)<0.4) rB1=1;
    else if (deltaR(etaB1, etaB4_gen, phiB1, phiB4_gen)<0.4) rB1=1;

    if (deltaR(etaB2, etaB1_gen, phiB2, phiB1_gen)<0.4) rB2=1;
    else if (deltaR(etaB2, etaB2_gen, phiB2, phiB2_gen)<0.4) rB2=1;
    else if (deltaR(etaB2, etaB3_gen, phiB2, phiB3_gen)<0.4) rB2=1;
    else if (deltaR(etaB2, etaB4_gen, phiB2, phiB4_gen)<0.4) rB2=1;

    outtree->Fill();
  }
  
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
