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
//#include "HttStyles.h"
#endif

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;

Float_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 );

//void wp(TString infilename = "hh.root", Float_t etaLo=0, Float_t etaHi=1.5, Float_t ptLo=30, Float_t ptHi=900, Float_t tarEff=0.97)
void wp(TString infilename = "hh.root", Float_t etaLo=1.5, Float_t etaHi=4.0, Float_t ptLo=30, Float_t ptHi=900, Float_t tarEff=0.97)
{

  TH1D* vBetaStarHard;
  TH1D* vMeanSqDeltaRHard;
  TH1D* vBetaStarPU;
  TH1D* vMeanSqDeltaRPU;

  char histname[50];

  sprintf(histname, "betaStarHard%i%i", etaLo, ptLo);
  vBetaStarHard = new TH1D(histname, histname, 100, 0, 1);
  sprintf(histname, "betaStarPU%i%i", etaLo, ptLo);
  vBetaStarPU = new TH1D(histname, histname, 100, 0, 1);
  sprintf(histname, "meanSqDeltaRHard%i%i", etaLo, ptLo);
  vMeanSqDeltaRHard = new TH1D(histname, histname, 25, 0, 0.25);
  sprintf(histname, "meanSqDeltaRPU%i%i", etaLo, ptLo);
  vMeanSqDeltaRPU = new TH1D(histname, histname, 25, 0, 0.25);
  
  Float_t pt=0, eta=0, phi=0;
  Float_t beta=0, betastar=0, mdrsq=0;
  UInt_t pu=0;
  
  TFile *infile = new TFile(infilename, "READ");
  TTree *intree = (TTree*) infile->Get("Events");

  intree->SetBranchAddress("pt", &pt);
  intree->SetBranchAddress("eta", &eta);
  intree->SetBranchAddress("phi", &phi);
  intree->SetBranchAddress("beta", &beta);
  intree->SetBranchAddress("betastar", &betastar);
  intree->SetBranchAddress("mdrsq", &mdrsq);
  intree->SetBranchAddress("pu", &pu);

  Float_t totalHard=0;
  Float_t totalPU=0;
  Float_t cutHard=0;
  Float_t cutPU=0;

  for (Int_t iEntry=0; iEntry<intree->GetEntries(); iEntry++) { // entry loop
    intree->GetEntry(iEntry);

    if ( (fabs(eta)<etaLo) || (fabs(eta)>etaHi) ) continue;
    if ( (pt<ptLo) || (pt>ptHi) ) continue;

    if (pu==1) {
      totalPU++;
      vBetaStarPU->Fill(betastar);
      vMeanSqDeltaRPU->Fill(mdrsq);
    }
    else {
      totalHard++;
      vBetaStarHard->Fill(betastar);
      vMeanSqDeltaRHard->Fill(mdrsq);

    }

  } // end event loop

  Float_t eff=0;
  Float_t rej=0;
  
  Float_t best_eff=0;
  Float_t best_rej=1;
  Float_t best_betas=0;

  //cout << "betaStar " << vBetaStarHard->GetEntries() << ", " << vBetaStarPU->GetEntries() << endl;

  for (Int_t i=0; i<vBetaStarHard->GetNbinsX()+2; i++) {
    eff = vBetaStarHard->Integral(0,i)/vBetaStarHard->GetEntries();
    rej = vBetaStarPU->Integral(0,i)/vBetaStarPU->GetEntries();

    //if ((eff>0.9)&&(rej<1.0)) cout << vBetaStarHard->GetBinLowEdge(i)+vBetaStarHard->GetBinWidth(i) << ", " << rej << ", " << eff << endl;

    if ((eff>tarEff)&&(eff>best_eff)&&(rej<best_rej)) {
      best_eff=eff;
      best_rej=rej;
      best_betas=vBetaStarHard->GetBinLowEdge(i)+vBetaStarHard->GetBinWidth(i);
    }
    
  }
  //cout << endl;
  //cout << best_eff << ", " << best_rej << ", " << best_betas << endl;

  //cout << "meanSqDeltaR " << vMeanSqDeltaRHard->GetEntries() << ", " << vMeanSqDeltaRPU->GetEntries() << endl;

  best_eff=0;
  best_rej=1;
  Float_t best_r=0;

  for (Int_t i=0; i<vMeanSqDeltaRHard->GetNbinsX()+1; i++) {
    eff = vMeanSqDeltaRHard->Integral(0,i)/vMeanSqDeltaRHard->GetEntries();
    rej = vMeanSqDeltaRPU->Integral(0,i)/vMeanSqDeltaRPU->GetEntries();
    
    //if ((eff>0.9)&&(rej<1.0)) cout << vMeanSqDeltaRHard->GetBinLowEdge(i)+vMeanSqDeltaRHard->GetBinWidth(i) << ", " << rej << ", " << eff << endl;

    if ((eff>tarEff)&&(eff>best_eff)&&(rej<best_rej)) {
      best_eff=eff;
      best_rej=rej;
      best_r=vMeanSqDeltaRHard->GetBinLowEdge(i)+vMeanSqDeltaRHard->GetBinWidth(i);
    }
  }
  
  //cout << endl;
  //cout << best_eff << ", " << best_rej << ", " << best_r << endl;

  totalPU=0;
  totalHard=0;

  for (Int_t iEntry=0; iEntry<intree->GetEntries(); iEntry++) { // entry loop
    intree->GetEntry(iEntry);

    if ( (fabs(eta)<etaLo) || (fabs(eta)>etaHi) ) continue;
    if ( (pt<ptLo) || (pt>ptHi) ) continue;

    if (pu==1) {
      totalPU++;
      if ((betastar<best_betas)&& (mdrsq<best_r)) cutPU++;
    }
    else {
      totalHard++;
      if ((betastar<best_betas)&& (mdrsq<best_r)) cutHard++;
    }

  } // end event loop

  cout << endl;
  cout << etaLo << " < |eta| < " << etaHi << endl;
  cout << ptLo << " < pt < " << ptHi << endl;
  cout << "Total events " << totalHard << "+" << totalPU << "=" << totalHard+totalPU << endl;
  cout << "BetaS:    " << best_betas << endl;
  cout << "MeanSqDR: " << best_r << endl;
  cout << "Eff:      " << cutHard/totalHard << " +- " << cutHard/totalHard*TMath::Sqrt(1/cutHard+1/totalHard) << endl;
  cout << "Rej:      " << cutPU/totalPU << " +- " << cutPU/totalPU*TMath::Sqrt(1/cutPU+1/totalPU) << endl;
  
}

Float_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 ) {

  const Float_t pi = 3.14159265358979;

  Float_t etaDiff = (eta1-eta2);
  Float_t phiDiff = fabs(phi1-phi2);
  while ( phiDiff>pi ) phiDiff = fabs(phiDiff-2.0*pi);

  Float_t deltaRSquared = etaDiff*etaDiff + phiDiff*phiDiff;

  return TMath::Sqrt(deltaRSquared);

}
