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
#include <THStack.h>
#include <TLine.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <vector>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <TGaxis.h>
#include <TLorentzVector.h>
#include "Math/LorentzVector.h"
#include "yieldData.hh"
#endif

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;

void confParse(const TString conf, vector<TString> &sampleNames, vector<TString> &sampleTitles, vector<Int_t> &sampleColors);

Float_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 );
Double_t deltaPhi(const Double_t phi1, const Double_t phi2);

yieldData getYield(const TString infilename, Float_t TAU_PT_MIN_H, Float_t TAU_PT_MIN_L, Float_t TAU_PT_MIN_HH, Float_t JET_PT_MIN, Float_t TAU_ETA_MAX, Float_t JET_ETA_MAX);

void vbf() {

  Float_t TAU_PT_MIN_H = 30;
  Float_t TAU_PT_MIN_L = 30;
  Float_t TAU_PT_MIN_HH = 30;
  Float_t JET_PT_MIN = 30;
  
  Float_t TAU_ETA_MAX = 4.0;
  Float_t JET_ETA_MAX = 4.7;

  yieldData signal = getYield("/afs/cern.ch/work/j/jlawhorn/public/vbf/ntuples/Bjj-vbf.root", TAU_PT_MIN_H, TAU_PT_MIN_L, TAU_PT_MIN_HH, JET_PT_MIN, TAU_ETA_MAX, JET_ETA_MAX);
  yieldData bkgd = getYield("/afs/cern.ch/work/j/jlawhorn/public/vbfntuples/vbf_bgd.root", TAU_PT_MIN_H, TAU_PT_MIN_L, TAU_PT_MIN_HH, JET_PT_MIN, TAU_ETA_MAX, JET_ETA_MAX);

  cout << setprecision(1) << fixed;

  cout << "SIGNAL " << setw(10) << "tt" << setw(10) << "mt" << setw(10)  << "et" << setw(10) << "em" << endl;
  cout << "  All  " << setw(10) << signal.total_tt*3e3 << setw(10) << signal.total_mt*3e3 << setw(10) << signal.total_et*3e3 << setw(10) << signal.total_em*3e3 << endl;
  cout << " Tight " << setw(10) << "--" << setw(10) << signal.mt_tight*3e3 << setw(10) << signal.et_tight*3e3 << setw(10) << signal.em_tight*3e3 << endl;
  cout << " Loose " << setw(10) << "--" << setw(10) << signal.mt_loose*3e3 << setw(10) << signal.et_loose*3e3 << setw(10) << signal.em_loose*3e3 << endl;

  cout << " BKGD  " << setw(10) << "tt" << setw(10) << "mt" << setw(10)  << "et" << setw(10) << "em" << endl;
  cout << "  All  " << setw(10) << bkgd.total_tt*3e3 << setw(10) << bkgd.total_mt*3e3 << setw(10) << bkgd.total_et*3e3 << setw(10) << bkgd.total_em*3e3 << endl;
  cout << " Tight " << setw(10) << "--" << setw(10) << bkgd.mt_tight*3e3 << setw(10) << bkgd.et_tight*3e3 << setw(10) << bkgd.em_tight*3e3 << endl;
  cout << " Loose " << setw(10) << "--" << setw(10) << bkgd.mt_loose*3e3 << setw(10) << bkgd.et_loose*3e3 << setw(10) << bkgd.em_loose*3e3 << endl;

}

yieldData getYield(const TString infilename, Float_t TAU_PT_MIN_H, Float_t TAU_PT_MIN_L, Float_t TAU_PT_MIN_HH, Float_t JET_PT_MIN, Float_t TAU_ETA_MAX, Float_t JET_ETA_MAX) {

  // tau decay modes
  enum { hadron=1, electron, muon };

  // tautau decay modes
  enum { all=0, dijet, jetmu, jetele, muele };

  yieldData data;

  data.total_yield=0;
  data.total_tt=0; data.total_mt=0; data.total_et=0; data.total_em=0;
  data.mt_tight=0; data.et_tight=0; data.em_tight=0;
  data.mt_loose=0; data.et_loose=0; data.em_loose=0;

  Float_t eventWeight=1;
  UInt_t eventType;
  UInt_t tauCat1, tauCat2;
  UInt_t tFake1, tFake2;
  Float_t met, metPhi;
  Float_t dEta, mJJ;

  LorentzVector *sRecoJet1=0, *sRecoJet2=0;
  LorentzVector *sRecoTau1=0, *sRecoTau2=0;
  LorentzVector *sGenTau1=0, *sGenTau2=0;
  LorentzVector *sGenJetTau1=0, *sGenJetTau2=0;

  TFile *infile;
  TTree *intree;

  //cout << "Processing  " << infilename << " ..." << endl;
  infile = new TFile(infilename); assert(infile);
  intree = (TTree*) infile->Get("Events"); assert(intree);
  
  intree->SetBranchAddress("eventWeight",    &eventWeight);
  intree->SetBranchAddress("eventType",      &eventType);
  intree->SetBranchAddress("tauCat1",        &tauCat1);
  intree->SetBranchAddress("tauCat2",        &tauCat2);
  intree->SetBranchAddress("tFake1",         &tFake1);
  intree->SetBranchAddress("tFake2",         &tFake2);
  intree->SetBranchAddress("met",            &met);
  intree->SetBranchAddress("metPhi",         &metPhi);
  intree->SetBranchAddress("dEta",           &dEta);
  intree->SetBranchAddress("mJJ",            &mJJ);
  intree->SetBranchAddress("sGenTau1",       &sGenTau1);      // 4-vector for generator leading tau
  intree->SetBranchAddress("sGenTau2",       &sGenTau2);      // 4-vector for generator second tau
  intree->SetBranchAddress("sGenJetTau1",    &sGenJetTau1);   // 4-vector for generator leading tau
  intree->SetBranchAddress("sGenJetTau2",    &sGenJetTau2);   // 4-vector for generator second tau
  intree->SetBranchAddress("sRecoTau1",      &sRecoTau1);     // 4-vector for reconstructed leading tau
  intree->SetBranchAddress("sRecoTau2",      &sRecoTau2);     // 4-vector for reconstructed second tau
  intree->SetBranchAddress("sRecoJet1",      &sRecoJet1);       // 4-vector for reconstructed leading b-jet
  intree->SetBranchAddress("sRecoJet2",      &sRecoJet2);       // 4-vector for reconstructed second b-jet

  Int_t CAT=-1, BOOST=-1;
  
  for(UInt_t iEntry=0; iEntry<intree->GetEntries(); iEntry++) { // entry loop
    intree->GetEntry(iEntry);

    //eventWeight=1;

    CAT=-1; BOOST=-1;

    // skip events that don't have 2 reco jet's and 2 reco tau's
    if ( (sRecoJet1->Pt()==999) || (sRecoJet2->Pt()==999) ) continue;
    if ( (sRecoTau1->Pt()==999) || (sRecoTau2->Pt()==999) ) continue;
    
    data.total_yield+=eventWeight;

    if ((sRecoJet1->Pt()<JET_PT_MIN)||(sRecoJet2->Pt()<JET_PT_MIN)) continue;
    if ((fabs(sRecoJet1->Eta())>JET_ETA_MAX)||(fabs(sRecoJet2->Eta())>JET_ETA_MAX)) continue;

    if ((fabs(sRecoTau1->Eta())>TAU_ETA_MAX)||(fabs(sRecoTau2->Eta())>TAU_ETA_MAX)) continue;

    if ( ( (tauCat1==electron) && (tauCat2==electron) ) || ( (tauCat1==muon) && (tauCat2==muon) ) ) continue;
    else if ( ( tauCat1 == hadron ) && ( tauCat2==hadron ) ) {
      if ( (sRecoTau1->Pt()<TAU_PT_MIN_H)||(sRecoTau2->Pt()<TAU_PT_MIN_HH) ) continue;
      CAT=2;
    }
    else if ( ( tauCat1==hadron ) && ( tauCat2==muon) ) {
      if ( (sRecoTau1->Pt()<TAU_PT_MIN_H)||(sRecoTau2->Pt()<TAU_PT_MIN_L) ) continue;
      CAT=3;
    }
    else if ( ( tauCat2==hadron ) && ( tauCat1==muon) ) {
      if ( (sRecoTau1->Pt()<TAU_PT_MIN_L)||(sRecoTau2->Pt()<TAU_PT_MIN_H) ) continue;
      CAT=3;
    }
    else if ( (tauCat1==hadron) && (tauCat2==electron) ) {
      if ( (sRecoTau1->Pt()<TAU_PT_MIN_H)||(sRecoTau2->Pt()<TAU_PT_MIN_L) ) continue;
      CAT=4;
    }
    else if ( (tauCat2==hadron) && (tauCat1==electron) ) {
      if ( (sRecoTau1->Pt()<TAU_PT_MIN_L)||(sRecoTau2->Pt()<TAU_PT_MIN_H) ) continue;
      CAT=4;
    }
    else if ( ( (tauCat1==muon) && (tauCat2==electron) ) || ( (tauCat2==muon) && (tauCat1==electron) ) ) {
      if ( (sRecoTau1->Pt()<TAU_PT_MIN_L)||(sRecoTau2->Pt()<TAU_PT_MIN_L) ) continue;
      CAT=5;
    }

    LorentzVector hTT=*sRecoTau1+*sRecoTau2;

    if ((hTT.M()<50)||(hTT.M()>150)) continue;

    if (mJJ<500) continue;

    if (fabs(dEta)<3.5) continue;


    // TAU-TAU
    if (CAT==2) {
      data.total_tt+=eventWeight; 
    }

    // MU-TAU
    else if (CAT==3) {
      data.total_mt+=eventWeight; 
      if ((mJJ>700)&&(fabs(dEta)>4.0)) data.mt_tight+=eventWeight;
      else data.mt_loose+=eventWeight;
    }

    // ELE-TAU
    else if (CAT==4) {
      data.total_et+=eventWeight; 
      if ((mJJ>700)&&(fabs(dEta)>4.0)) data.et_tight+=eventWeight;
      else data.et_loose+=eventWeight;
    }

    // ELE-MU
    else if (CAT==5) {
      data.total_em+=eventWeight; 
      if ((mJJ>700)&&(fabs(dEta)>4.0)) data.em_tight+=eventWeight;
      else data.em_loose+=eventWeight;
    }

  } // end entry loop
  
  delete infile;
  infile=0, intree=0;

  return data;

}

void confParse(const TString conf, 
	       vector<TString> &sampleNames, 
	       vector<TString> &sampleTitles, 
	       vector<Int_t> &sampleColors) {

  ifstream ifs;
  ifs.open(conf.Data()); assert(ifs.is_open());
  string line;
  while(getline(ifs,line)) {

    if( (line[0]=='#') || (line[0]==' ') ) continue;

    string fname;
    string title;
    Int_t color;
    stringstream ss(line);
    ss >> fname >> title >> color;
    sampleNames.push_back(fname);
    sampleTitles.push_back(title);
    sampleColors.push_back(color);

  }
  ifs.close();

}

Double_t deltaPhi(Double_t phi1, Double_t phi2) 
{
  // Compute dPhi between two given angles. Results is in [0,pi].
  const Double_t pi = 3.14159265358979;
  Double_t dphi = fabs(phi1-phi2);
  while (dphi>pi)
    dphi = fabs(dphi - 2.0*pi);

  return dphi;
}

Float_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 ) {

  const Float_t pi = 3.14159265358979;

  Float_t etaDiff = (eta1-eta2);
  Float_t phiDiff = fabs(phi1-phi2);
  while ( phiDiff>pi ) phiDiff = fabs(phiDiff-2.0*pi);

  Float_t deltaRSquared = etaDiff*etaDiff + phiDiff*phiDiff;

  return TMath::Sqrt(deltaRSquared);

}
