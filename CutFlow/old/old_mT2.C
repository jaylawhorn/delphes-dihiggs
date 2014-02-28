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
#include <TLine.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <vector>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include "Math/LorentzVector.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#endif

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;

Float_t eT(LorentzVector *particle);

TVector2 sumPT(LorentzVector *tau1, LorentzVector *tau2, Float_t met, Float_t metPhi);

Float_t transMass(LorentzVector *a, LorentzVector *b);

Float_t stransverseMass(LorentzVector *b1, LorentzVector *c1, LorentzVector *b2, LorentzVector *c2);

void mT2() {

  ROOT::Math::Minimizer *min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
  min->SetMaxFunctionCalls(1000000);
  min->SetTolerance(0.01);
  min->SetPrintLevel(1);
  ROOT::Math::Functor f(&toMin,2);

  Double_t step[2] = {0.01, 0.01};
  Double_t variable[2] = {-1, 1.2};
  
  min->SetFunction(f);
  min->SetVariable(0,"x", variable[0], step[0]);
  min->SetVariable(1,"y", variable[1], step[1]);

  min->Minimize();

  const Double_t *xs = min->X();
  cout << "Min: " << xs[0] << ", " << xs[1] << ", " << min->MinValue() << endl;

  /*
  Float_t eventWeight=1;
  Float_t met, metPhi;
  UInt_t bTag1, bTag2;
  UInt_t tauCat1, tauCat2;
  LorentzVector *recoB1=0, *recoB2=0;
  LorentzVector *recoTau1=0, *recoTau2=0;
  LorentzVector *recoLeadJet=0, *recoExtraJet=0;
  LorentzVector *nu1=0, *nu2=0;

  Float_t nu1Pt=0, nu1Eta=0, nu1Phi=0, nu1M=0;
  Float_t nu2Pt=0, nu2Eta=0, nu2Phi=0, nu2M=0;

  ROOT::Math::Minimizer *min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
  min->SetMaxFunctionCalls(1000000);
  min->SetTolerance(0.01);
  min->SetPrintLevel(1);

  ROOT::Math::Functor f(&stransverseMass, 4);

  TFile *infile;
  TTree *intree;

  TH1D *hist = new TH1D("hist", "hist",50,0,300);

  TString infilename = "/afs/cern.ch/work/k/klawhorn/delphes_working/HHToTTBB_14TeV.root";
  infile = new TFile(infilename); assert(infile);
  intree = (TTree*) infile->Get("Events"); assert(intree);

  intree->SetBranchAddress("eventWeight",    &eventWeight);
  intree->SetBranchAddress("tauCat1",        &tauCat1);
  intree->SetBranchAddress("tauCat2",        &tauCat2);
  intree->SetBranchAddress("bTag1",          &bTag1);
  intree->SetBranchAddress("bTag2",          &bTag2);
  intree->SetBranchAddress("met",            &met);
  intree->SetBranchAddress("metPhi",         &metPhi);
  intree->SetBranchAddress("recoTau1",       &recoTau1);     // 4-vector for reconstructed leading tau
  intree->SetBranchAddress("recoTau2",       &recoTau2);     // 4-vector for reconstructed second tau
  intree->SetBranchAddress("recoB1",         &recoB1);       // 4-vector for reconstructed leading b-jet
  intree->SetBranchAddress("recoB2",         &recoB2);       // 4-vector for reconstructed second b-jet
  intree->SetBranchAddress("recoLeadJet",    &recoLeadJet);  // 4-vector for reconstructed leading jet
  intree->SetBranchAddress("recoExtraJet",   &recoExtraJet); // 4-vector for reconstructed extra jet
  
  for(UInt_t iEntry=0; iEntry<intree->GetEntries(); iEntry++) { // entry loop
    intree->GetEntry(iEntry);

    TVector2 ptConstraint = sumPT(recoTau1, recoTau2, met, metPhi);

  } // end entry loop

  leadingB->GetYaxis()->SetRangeUser(0,250);
  leadingB->Draw();
  secondB->SetLineColor(kRed);
  secondB->Draw("same");
  leadingT->SetLineColor(kGreen);
  leadingT->Draw("same");
  secondT->SetLineColor(kBlue);
  secondT->Draw("same");

  delete infile;
  infile=0, intree=0;
  */
}

Float_t eT(LorentzVector *particle) {

  return TMath::Sqrt(particle->M2() + particle->Perp2());

}

TVector2 sumPT( LorentzVector *tau1, LorentzVector *tau2, Float_t met, Float_t metPhi) {

  TVector2 tauT1(0,0);
  TVector2 tauT2(0,0);
  TVector2 metT(0,0);

  tauT1.SetMagPhi(tau1->Pt(), tau1->Phi());
  tauT2.SetMagPhi(tau2->Pt(), tau2->Phi());
  metT.SetMagPhi(met, metPhi);

  TVector2 sumPT = tauT1+tauT2+metT;
  return sumPT;

}

Float_t transMass(LorentzVector *a, LorentzVector *b) {

  Float_t mass2 = a->M2()+b->M2()+2*(eT(a)*eT(b)-a->Dot(*b));

  return TMath::Sqrt(mass2);

}

Float_t stransverseMass(LorentzVector *b1, LorentzVector *c1, LorentzVector *b2, LorentzVector *c2) {

  Float_t mT1 = transMass(b1, c1);
  Float_t mT2 = transMass(b2, c2);

  return TMath::Max(mT1, mT2);

}
