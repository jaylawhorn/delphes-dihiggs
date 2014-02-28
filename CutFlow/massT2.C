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
#include "../MitStyleRemix.hh"
#include "../CPlot.hh"
#endif

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;

Double_t mTsq(TVector2 bT, TVector2 cT, Float_t mB, Float_t mC);

class smT2
{
public:
  smT2() { name="thing"; }
  double operator()(const double *thing);

  void SetB1(TVector2 b1) { B1=b1; }
  void SetB2(TVector2 b2) { B2=b2; }
  void SetT1(TVector2 t1) { T1=t1; }
  void SetT2(TVector2 t2) { T2=t2; }
  void SetMPT(TVector2 mpt) { MPT=mpt; }
  void SetMB1(Float_t mB1) { MB1=mB1; }
  void SetMB2(Float_t mB2) { MB2=mB2; }
  void SetMT1(Float_t mT1) { MT1=mT1; }
  void SetMT2(Float_t mT2) { MT2=mT2; }

protected:
  TString name;
  TVector2 B1;
  TVector2 B2;
  TVector2 T1;
  TVector2 T2;
  TVector2 MPT;
  Float_t MB1;
  Float_t MB2;
  Float_t MT1;
  Float_t MT2;
};

double smT2::operator()(const double *thing) {
  double cT = thing[0];
  double cPhi = thing[1];
  TVector2 temp1(0,0);
  temp1.SetMagPhi(cT,cPhi);
  TVector2 temp2 = MPT-temp1;
  SetT1(temp1);
  SetT2(temp2);

  return TMath::Max(TMath::Sqrt(mTsq(B1, T1, MB1, MT1)), TMath::Sqrt(mTsq(B2, T2, MB2, MT2)));
}

void massT2() {

  ROOT::Math::Minimizer *min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
  // set tolerance , etc...
  min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2 
  min->SetTolerance(10.0);
  min->SetPrintLevel(0);

  Float_t eventWeight=1;
  Float_t met, metPhi;
  LorentzVector *recoB1=0, *recoB2=0;
  LorentzVector *recoTau1=0, *recoTau2=0;

  TVector2 tau1(0,0), tau2(0,0), mpt(0,0);
  TVector2 b1(0,0), b2(0,0);
  Float_t mTau1=0, mTau2=0;
  Float_t mB1=0, mB2=0;

  Double_t mt2=0;

  TH1D *hist = new TH1D("hist", "hist", 50, 0, 600); hist->Sumw2();

  //TString infilename = "/afs/cern.ch/work/k/klawhorn/delphes_working/HHToTTBB_14TeV.root";
  TString infilename = "/afs/cern.ch/work/k/klawhorn/delphes_working/tt_14TEV.root";
  TFile *infile = new TFile(infilename); assert(infile);
  TTree *intree = (TTree*) infile->Get("Events"); assert(intree);

  intree->SetBranchAddress("eventWeight",    &eventWeight);
  intree->SetBranchAddress("met",            &met);
  intree->SetBranchAddress("metPhi",         &metPhi);
  intree->SetBranchAddress("recoTau1",       &recoTau1);     // 4-vector for reconstructed leading tau
  intree->SetBranchAddress("recoTau2",       &recoTau2);     // 4-vector for reconstructed second tau
  intree->SetBranchAddress("recoB1",         &recoB1);       // 4-vector for reconstructed leading b-jet
  intree->SetBranchAddress("recoB2",         &recoB2);       // 4-vector for reconstructed second b-jet

  for(UInt_t iEntry=0; iEntry<200000/*intree->GetEntries()*/; iEntry++) { // entry loop
    intree->GetEntry(iEntry);

    tau1.SetMagPhi(recoTau1->Pt(), recoTau1->Phi());
    tau2.SetMagPhi(recoTau2->Pt(), recoTau2->Phi());
    mTau1=recoTau1->M();
    mTau2=recoTau2->M();

    b1.SetMagPhi(recoB1->Pt(), recoB1->Phi());
    b2.SetMagPhi(recoB2->Pt(), recoB2->Phi());
    mB1=recoB1->M();
    mB2=recoB2->M();

    mpt.SetMagPhi(met, metPhi);

    TVector2 sumPt = tau1+tau2+mpt;

    smT2 testing = smT2();
    testing.SetB1(b1);
    testing.SetB2(b2);
    testing.SetMPT(sumPt);
    testing.SetMB1(mB1);
    testing.SetMB2(mB2);
    testing.SetMT1(mTau1);
    testing.SetMT2(mTau2);

    TVector2 c1=sumPt;
    TVector2 c2=sumPt-c1;

    //testing(c1.Mod(), c1.Phi());
    ROOT::Math::Functor f(testing,2);
    double step[2] = {0.1, 0.1};
    double variable[2] = { 0.5*c1.Mod(), 0.0 };

    //cout << variable[0] << " " << variable[1] << " " <<testing(variable) << endl;

    min->SetFunction(f);
    min->SetLimitedVariable(0,"cT",variable[0], step[0], 0.0, sumPt.Mod());
    min->SetLimitedVariable(1,"cPhi",variable[1], step[1], 0.0, TMath::Pi());

    min->Minimize();
    const double *xs = min->X();
    std::cout << "Minimum: f(" << xs[0] << "," << xs[1] << "): " << min->MinValue()  << std::endl;

    mt2 = min->MinValue();

    hist->Fill(mt2, eventWeight);

  }

  TCanvas *c1 = MakeCanvas("c1", "c1", 800, 600);

  Float_t scale = hist->Integral();
  hist->Scale(1/scale);

  hist->SetTitle("t#bar{t}");
  hist->GetXaxis()->SetTitle("mT2");
  hist->GetYaxis()->SetTitle("Normalized Events");

  hist->Draw("hist");

  c1->SaveAs("tt_mt2.png");

  delete infile;
  infile=0, intree=0;

}

Double_t mTsq(TVector2 bT, TVector2 cT, Float_t mB, Float_t mC) {
  Float_t eB = TMath::Sqrt(mB*mB+ bT*bT);
  Float_t eC = TMath::Sqrt(mC*mC+ cT*cT);

  return mB*mB+mC*mC+2*(eB*eC - bT*cT);
}
