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
#include <TLegend.h>
#include <vector>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include "Math/LorentzVector.h"

#include "modules/Delphes.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "classes/DelphesClasses.h"

#include "MitStyleRemix.hh"
#endif

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;

Float_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 );

//void electrons(const TString inputfile="root://eoscms.cern.ch//store/group/upgrade/delphes/ProdJun14/tt-4p-0-600-v1510_14TEV/tt-4p-0-600-v1510_14TEV_216084397_PhaseII_Conf4_140PileUp_seed216084402_5of5.root") {
//void electrons(const TString inputfile="root://eoscms.cern.ch//store/group/upgrade/delphes/ProdJun14/tt-4p-1100-1700-v1510_14TEV/tt-4p-1100-1700-v1510_14TEV_204948586_PhaseII_Conf4_140PileUp_seed204948588_2of5.root") {
void electrons(const TString input="hgg.txt") {

  // read input input file                                                                                              
  TChain chain("Delphes");
  ifstream ifs;
  ifs.open(input.Data());
  assert(ifs.is_open());
  string line;
  getline(ifs,line);
  while(getline(ifs,line)) { 
    chain.Add(TString(line));
  }
  ifs.close();
  
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  TClonesArray *branchPhoton = treeReader->UseBranch("Photon");  
  //TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  //TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  
  Int_t kMatch=0;
  //Electron *ele;
  //Electron *matched;
  Photon *pho;
  //Photon *matched;
  //Muon *muon;
  //Muon *matched;
  GenParticle *part;
  
  /*TProfile *hPt = new TProfile("hPt", "hPt", 12, 0, 300);
  TProfile *hEta = new TProfile("hEta", "hEta", 16, -4.0, 4.0);
  TProfile *hPtIso4 = new TProfile("hPtIso4", "hPtIso4", 12, 0, 300);
  TProfile *hEtaIso4 = new TProfile("hEtaIso4", "hEtaIso4", 16, -4.0, 4.0);
  TProfile *hPtIso2 = new TProfile("hPtIso2", "hPtIso2", 12, 0, 300);
  TProfile *hEtaIso2 = new TProfile("hEtaIso2", "hEtaIso2", 16, -4.0, 4.0);
  TProfile *hPtIso0 = new TProfile("hPtIso0", "hPtIso0", 12, 0, 300);
  TProfile *hEtaIso0 = new TProfile("hEtaIso0", "hEtaIso0", 16, -4.0, 4.0);
  TProfile *hPtIso2m = new TProfile("hPtIso2m", "hPtIso2m", 12, 0, 300);
  TProfile *hEtaIso2m = new TProfile("hEtaIso2m", "hEtaIso2m", 16, -4.0, 4.0);
  TProfile *hPtIso4m = new TProfile("hPtIso4m", "hPtIso4m", 12, 0, 300);
  TProfile *hEtaIso4m = new TProfile("hEtaIso4m", "hEtaIso4m", 16, -4.0, 4.0);*/

  //TH1D *ptdist = new TH1D("ptdist", "ptdist", 10, 0, 100);

  TH1D *hIso = new TH1D("hIso", "hIso", 16, -1, 3);
  
  for (Int_t i=0; i<numberOfEntries; i++) {
    treeReader->ReadEntry(i);

    for (Int_t j=0; j<branchPhoton->GetEntries(); j++) {
      pho = (Photon*) branchPhoton->At(j);

      hIso->Fill(pho->IsolationVar);
      
    }
    
    /*    for (Int_t j=0; j<branchParticle->GetEntries(); j++) {
      part = (GenParticle*) branchParticle->At(j);
      
      if (abs(part->PID) != 11) continue;
      //if (abs(part->PID) != 13) continue;
      matched=0;
      Int_t id=0;
      for (Int_t k=0; k<branchElectron->GetEntries(); k++) {
	ele = (Electron*) branchElectron->At(k);
	if (deltaR(ele->Eta, part->Eta, ele->Phi, part->Phi) < 0.4) {
	  matched = (Electron*) branchElectron->At(kMatch);
	  id=1;
	}
	}*/
      /*for (Int_t k=0; k<branchMuon->GetEntries(); k++) {
	muon = (Muon*) branchMuon->At(k);
	if (deltaR(muon->Eta, part->Eta, muon->Phi, part->Phi) < 0.4) {
	  matched = (Muon*) branchMuon->At(kMatch);
	  id=1;
	}
	}*/
      
    /*      if (matched && part->PT>30) {
      //if (matched) {
	//if (matched->IsolationVar>0.4 && part->Eta>3.0) {
	if (matched->IsolationVar>0.4) {
	  hPtIso4->Fill(part->PT, 0);
	  hEtaIso4->Fill(part->Eta, 0);
	  hPtIso2->Fill(part->PT, 0);
	  hEtaIso2->Fill(part->Eta, 0);
	  hPtIso0->Fill(part->PT, 0);
	  hEtaIso0->Fill(part->Eta, 0);
	  hPtIso2m->Fill(part->PT, 0);
	  hEtaIso2m->Fill(part->Eta, 0);
	  hPtIso4m->Fill(part->PT, 0);
	  hEtaIso4m->Fill(part->Eta, 0);
	}
	//else if (matched->IsolationVar>0.2 && part->Eta>3.0) {
	else if (matched->IsolationVar>0.2) {
	  hPtIso4->Fill(part->PT, 1);
	  hEtaIso4->Fill(part->Eta, 1);
	  hPtIso2->Fill(part->PT, 0);
	  hEtaIso2->Fill(part->Eta, 0);
	  hPtIso0->Fill(part->PT, 0);
	  hEtaIso0->Fill(part->Eta, 0);
	  hPtIso2m->Fill(part->PT, 0);
	  hEtaIso2m->Fill(part->Eta, 0);
	  hPtIso4m->Fill(part->PT, 0);
	  hEtaIso4m->Fill(part->Eta, 0);
	}
	//else if (matched->IsolationVar>0.0 && part->Eta>3.0) {
	else if (matched->IsolationVar>0.0) {
	  hPtIso4->Fill(part->PT, 1);
	  hEtaIso4->Fill(part->Eta, 1);
	  hPtIso2->Fill(part->PT, 1);
	  hEtaIso2->Fill(part->Eta, 1);
	  hPtIso0->Fill(part->PT, 0);
	  hEtaIso0->Fill(part->Eta, 0);
	  hPtIso2m->Fill(part->PT, 0);
	  hEtaIso2m->Fill(part->Eta, 0);
	  hPtIso4m->Fill(part->PT, 0);
	  hEtaIso4m->Fill(part->Eta, 0);
	}
	//else if (matched->IsolationVar>-0.25 && part->Eta>3.0) {
	else if (matched->IsolationVar>-0.25) {
	  hPtIso4->Fill(part->PT, 1);
	  hEtaIso4->Fill(part->Eta, 1);
	  hPtIso2->Fill(part->PT, 1);
	  hEtaIso2->Fill(part->Eta, 1);
	  hPtIso0->Fill(part->PT, 1);
	  hEtaIso0->Fill(part->Eta, 1);
	  hPtIso2m->Fill(part->PT, 0);
	  hEtaIso2m->Fill(part->Eta, 0);
	  hPtIso4m->Fill(part->PT, 0);
	  hEtaIso4m->Fill(part->Eta, 0);
	}
	//else if (matched->IsolationVar>-0.4 && part->Eta>3.0) {
	else if (matched->IsolationVar>-0.4) {
	  hPtIso4->Fill(part->PT, 1);
	  hEtaIso4->Fill(part->Eta, 1);
	  hPtIso2->Fill(part->PT, 1);
	  hEtaIso2->Fill(part->Eta, 1);
	  hPtIso0->Fill(part->PT, 1);
	  hEtaIso0->Fill(part->Eta, 1);
	  hPtIso2m->Fill(part->PT, 1);
	  hEtaIso2m->Fill(part->Eta, 1);
	  hPtIso4m->Fill(part->PT, 0);
	  hEtaIso4m->Fill(part->Eta, 0);
	}
	//else if (part->Eta>3.0) {
	else {
	  hPtIso4->Fill(part->PT, 1);
	  hEtaIso4->Fill(part->Eta, 1);
	  hPtIso2->Fill(part->PT, 1);
	  hEtaIso2->Fill(part->Eta, 1);
	  hPtIso0->Fill(part->PT, 1);
	  hEtaIso0->Fill(part->Eta, 1);
	  hPtIso2m->Fill(part->PT, 1);
	  hEtaIso2m->Fill(part->Eta, 1);
	  hPtIso4m->Fill(part->PT, 1);
	  hEtaIso4m->Fill(part->Eta, 1);
	}
      }
      }*/
  }
  
  TCanvas *c = MakeCanvas("c", "c", 600, 600);

  hIso->SetTitle("");
  hIso->GetXaxis()->SetTitle("Photon.IsolationVar");
  hIso->GetYaxis()->SetTitle("Events");
  hIso->Draw("");
  c->SaveAs("photonIsolation.png");
  
  /*
  hPtIso4m->SetLineColor(kBlue);
  hPtIso4m->SetMarkerColor(kBlue);
  hPtIso2m->SetLineColor(kGreen);
  hPtIso2m->SetMarkerColor(kGreen);
  hPtIso0->SetLineColor(kRed);
  hPtIso0->SetMarkerColor(kRed);
  hPtIso2->SetLineColor(kMagenta);
  hPtIso2->SetMarkerColor(kMagenta);
  hPtIso4->SetLineColor(kBlack);
  hPtIso4->SetMarkerColor(kBlack);
  hPtIso4->GetYaxis()->SetRangeUser(0,1.1);
  hPtIso4->Draw("");
  hPtIso2->Draw("same");
  hPtIso0->Draw("same");
  hPtIso2m->Draw("same");
  hPtIso4m->Draw("same");
  c->SaveAs("ptiso.png");
  
  hEtaIso4m->SetLineColor(kBlue);
  hEtaIso4m->SetMarkerColor(kBlue);
  hEtaIso2m->SetLineColor(kGreen);
  hEtaIso2m->SetMarkerColor(kGreen);
  hEtaIso0->SetLineColor(kRed);
  hEtaIso0->SetMarkerColor(kRed);
  hEtaIso2->SetLineColor(kMagenta);
  hEtaIso2->SetMarkerColor(kMagenta);
  hEtaIso4->SetLineColor(kBlack);
  hEtaIso4->SetMarkerColor(kBlack);
  hEtaIso4->GetYaxis()->SetRangeUser(0,1.1);
  hEtaIso4->Draw("");
  hEtaIso2->Draw("same");
  hEtaIso0->Draw("same");
  hEtaIso2m->Draw("same");
  hEtaIso4m->Draw("same");
  c->SaveAs("etaiso.png");
  */  
}

Float_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 ) {

  const Float_t pi = 3.14159265358979;

  Float_t etaDiff = (eta1-eta2);
  Float_t phiDiff = fabs(phi1-phi2);
  while ( phiDiff>pi ) phiDiff = fabs(phiDiff-2.0*pi);

  Float_t deltaRSquared = etaDiff*etaDiff + phiDiff*phiDiff;

  return TMath::Sqrt(deltaRSquared);

}
