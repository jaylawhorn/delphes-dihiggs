#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TStyle.h>
#include <TLegend.h>
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

void makeHist(TString file, TH1D* hist);

void debug() {

  TString file0="/afs/cern.ch/work/j/jlawhorn/public/delphesTest/test_more0.root";
  TString file1="/afs/cern.ch/work/j/jlawhorn/public/delphesTest/test_more1.root";
  TString file2="/afs/cern.ch/work/j/jlawhorn/public/delphesTest/test_more2.root";
  TString file3="/afs/cern.ch/work/j/jlawhorn/public/delphesTest/test_more3.root";
  TString file4="/afs/cern.ch/work/j/jlawhorn/public/delphesTest/test_more4.root";
  TString file5="/afs/cern.ch/work/j/jlawhorn/public/delphesTest/test_more5.root";
  TString file6="/afs/cern.ch/work/j/jlawhorn/public/delphesTest/test_more6.root";
  TString file7="/afs/cern.ch/work/j/jlawhorn/public/delphesTest/test_more7.root";
  TString file8="/afs/cern.ch/work/j/jlawhorn/public/delphesTest/test_more8.root";
  TString file9="/afs/cern.ch/work/j/jlawhorn/public/delphesTest/new_test.root";

  Int_t nbins=10;
  Int_t low=0;
  Int_t high=1.0;

  TH1D *hist0 = new TH1D("iso 0.01", "", nbins, low, high); hist0->SetLineColor(kRed);
  TH1D *hist1 = new TH1D("iso 0.4", "", nbins, low, high);
  TH1D *hist2 = new TH1D("iso 2.0", "", nbins, low, high);
  TH1D *hist3 = new TH1D("iso 3.0", "", nbins, low, high);
  TH1D *hist4 = new TH1D("iso 9999", "", nbins, low, high); hist4->SetLineColor(kBlack);
  TH1D *hist5 = new TH1D("iso 0.8", "", nbins, low, high); hist5->SetLineColor(kBlue);
  TH1D *hist6 = new TH1D("iso 1.0", "", nbins, low, high); hist6->SetLineColor(kGreen);
  TH1D *hist7 = new TH1D("iso 1.2", "", nbins, low, high); hist7->SetLineColor(kMagenta);
  TH1D *hist8 = new TH1D("iso 1.5", "", nbins, low, high); hist8->SetLineColor(42);
  TH1D *hist9 = new TH1D("new", "", nbins, low, high); hist9->SetLineColor(42);

  makeHist(file0, hist0);
  makeHist(file1, hist1);
  makeHist(file5, hist5);
  makeHist(file6, hist6);
  makeHist(file7, hist7);
  makeHist(file8, hist8);
  makeHist(file2, hist2);
  makeHist(file3, hist3);
  makeHist(file4, hist4);
  makeHist(file9, hist9);

  gStyle->SetOptStat(0);

  hist0->SetTitle("");
  hist0->GetXaxis()->SetTitle("min(#Delta R(jet,b-quark))");
  hist0->GetYaxis()->SetTitle("Normalized Events");
  hist0->Draw("");
  hist4->Draw("same");
  //hist5->Draw("same");
  //hist6->Draw("same");
  //hist7->Draw("same");
  //hist8->Draw("same");
  hist9->Draw("same");

  TLegend *leg = new TLegend(0.7, 0.6, 0.9, 0.9);
  leg->SetFillColor(0);
  leg->SetShadowColor(0);
  leg->AddEntry(hist0, "Iso=0.01", "l");
  //leg->AddEntry(hist5, "Iso=0.8", "l");
  //leg->AddEntry(hist6, "Iso=1.0", "l");
  //leg->AddEntry(hist7, "Iso=1.2", "l");
  //leg->AddEntry(hist8, "Iso=1.5", "l");
  leg->AddEntry(hist4, "Iso=9999", "l");
  leg->AddEntry(hist9, "No Iso", "l");
  leg->Draw("same");

}

Float_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 ) {

  const Float_t pi = 3.14159265358979;

  Float_t etaDiff = (eta1-eta2);
  Float_t phiDiff = fabs(phi1-phi2);
  while ( phiDiff>pi ) phiDiff = fabs(phiDiff-2.0*pi);

  Float_t deltaRSquared = etaDiff*etaDiff + phiDiff*phiDiff;

  return TMath::Sqrt(deltaRSquared);

}

void makeHist(TString file, TH1D* hist) {

  // read input file
  TChain chain("Delphes");
  chain.Add(file);
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");

  //set up loop variables
  GenParticle *part=0;
  Jet *jet=0;

  Float_t minDR=100;
  Float_t dr;

  Int_t ii=0;

  for (Int_t iEntry=0; iEntry<numberOfEntries; iEntry++) { // entry loop
    treeReader->ReadEntry(iEntry);

    for (Int_t iJet=0; iJet<branchJet->GetEntries(); iJet++) {
      jet = (Jet*) branchJet->At(iJet);

      minDR=100; 

      if (jet->BTag>0) ii++;

      for (Int_t iPart=0; iPart<branchParticle->GetEntries(); iPart++) {
	part = (GenParticle*) branchParticle->At(iPart);

	if ( fabs(part->PID) != 5 ) continue;

	dr=deltaR(jet->Eta, part->Eta, jet->Phi, part->Phi);
	if ( dr < minDR ) minDR=dr;
      }

      hist->Fill(minDR);
      
    }
    
  } // end entry loop

  cout << hist->GetName() << " - # b-jets: " << ii << " +- " << TMath::Sqrt(ii) << endl;

  Float_t n=hist->GetEntries();
  hist->Scale(1.0/n);

}
