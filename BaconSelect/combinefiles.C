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

#endif

Double_t btag(Double_t pt, Double_t eta);

void combinefiles(const TString input="hh_proc.txt",
		  const TString outputfile="hh_proc.root") {

  // read input file
  TChain *eventChain = new TChain("Events");
  TChain *countChain = new TChain("Count");

  ifstream ifs;
  ifs.open(input.Data());
  assert(ifs.is_open());
  string line;
  getline(ifs,line);
  while(getline(ifs,line)) { 
    eventChain->Add(TString(line));
    countChain->Add(TString(line));
  }
  ifs.close();

  Int_t iHmatch1=-1, iHmatch2=-1, iHmatch3=-1, iHmatch4=-1;
  Int_t eventType; Int_t genInfo; Float_t eventWeight; Int_t sampleNo;
  Int_t isBBTT; Int_t isBBGG; Int_t isBBBB; Int_t isVBFTT; Int_t isVBF4B;
  Float_t met, metPhi; Float_t ppMet, ppMetPhi; Float_t pileupMet, pileupMetPhi;
  Int_t nCentral=0, nBtag=0; Int_t centB=0; Int_t nJets=0;

  Int_t tauCat1=0, tauCat2=0;
  Int_t bTag1=0, bTag2=0, bTag3=0, bTag4=0;
  Int_t jbTag_tt1=0, jbTag_tt2=0, jbTag_6j1=0, jbTag_6j2=0;

  Double_t mt2=0; Double_t mt2puppi=0; Double_t mt2pileup=0;
  Double_t m_sv=0; Double_t m_svpileup=0; Double_t m_svpuppi=0;
  
  Float_t ptTau1, ptTau2, ptTrk1, ptTrk2, ptB1, ptB2, ptB3, ptB4, ptG1, ptG2, ptJet_tt1, ptJet_tt2, ptJet_6j1, ptJet_6j2, tauIso1, tauIso2;
  Float_t etaTau1, etaTau2, etaTrk1, etaTrk2, etaB1, etaB2, etaB3, etaB4, etaG1, etaG2, etaJet_tt1, etaJet_tt2, etaJet_6j1, etaJet_6j2;
  Float_t phiTau1, phiTau2, phiTrk1, phiTrk2, phiB1, phiB2, phiB3, phiB4, phiG1, phiG2, phiJet_tt1, phiJet_tt2, phiJet_6j1, phiJet_6j2;
  Float_t mTau1, mTau2, mTrk1, mTrk2, mB1, mB2, mB3, mB4, eG1, eG2, mJet_tt1, mJet_tt2, mJet_6j1, mJet_6j2;

  Float_t ptTau1_gen, ptTau2_gen, ptB1_gen, ptB2_gen, ptB3_gen, ptB4_gen, ptG1_gen, ptG2_gen, ptJet_tt1_gen, ptJet_tt2_gen, ptJet_6j1_gen, ptJet_6j2_gen;
  Float_t etaTau1_gen, etaTau2_gen, etaB1_gen, etaB2_gen, etaB3_gen, etaB4_gen, etaG1_gen, etaG2_gen, etaJet_tt1_gen, etaJet_tt2_gen, etaJet_6j1_gen, etaJet_6j2_gen;
  Float_t phiTau1_gen, phiTau2_gen, phiB1_gen, phiB2_gen, phiB3_gen, phiB4_gen, phiG1_gen, phiG2_gen, phiJet_tt1_gen, phiJet_tt2_gen, phiJet_6j1_gen, phiJet_6j2_gen;
  Float_t mTau1_gen, mTau2_gen, mB1_gen, mB2_gen, mB3_gen, mB4_gen, eG1_gen, eG2_gen, mJet_tt1_gen, mJet_tt2_gen, mJet_6j1_gen, mJet_6j2_gen;

  Float_t gamIso1, gamIso2;

  Float_t ptH1_gen, etaH1_gen, phiH1_gen, mH1_gen;
  Float_t ptH2_gen, etaH2_gen, phiH2_gen, mH2_gen;

  Float_t ptTau1_genJet, ptTau2_genJet, etaTau1_genJet, etaTau2_genJet, phiTau1_genJet, phiTau2_genJet, mTau1_genJet, mTau2_genJet;

  Float_t ptTT, ptBB1, ptBB2, ptGG, ptJJ_tt, ptJJ_6j, ptHH;
  Float_t etaTT, etaBB1, etaBB2, etaGG, etaJJ_tt, etaJJ_6j, etaHH;
  Float_t phiTT, phiBB1, phiBB2, phiGG, phiJJ_tt, phiJJ_6j, phiHH;
  Float_t mTT, mBB1, mBB2, mGG, mJJ_tt, mJJ_6j, mHH;

  Float_t dEta_tt=0, dEta_6j=0; Int_t n;

  Float_t rho_0=0, rho_1=0, rho_2=0;

  Float_t dEtaBB1=0, dEtaTT=0, dEtaHH=0;
  Float_t dPhiBB1=0, dPhiTT=0, dPhiHH=0;
  Float_t dRBB1=0, dRTT=0, dRHH=0;

  Float_t bdtVal=0;

  Float_t chi2=0, corrMet=0, corrMetPhi=0;
  Int_t conv=-999;

  TFile *outFile = new TFile(outputfile, "RECREATE");

  countChain->SetBranchAddress("n", &n);

  Double_t nevents=0;

  for (Int_t i=0; i<countChain->GetEntries(); i++) {
    countChain->GetEntry(i);
    nevents+=n;
  }
  cout << nevents << endl;

  eventChain->SetBranchAddress("mTT",            &mTT);
  eventChain->SetBranchAddress("eventWeight",    &eventWeight);

  eventChain->SetBranchAddress("ptB1",           &ptB1);
  eventChain->SetBranchAddress("ptB2",           &ptB2);
  eventChain->SetBranchAddress("etaB1",           &etaB1);
  eventChain->SetBranchAddress("etaB2",           &etaB2);

  eventChain->GetEntry(0);

  TTree *outTree=(TTree*)eventChain->GetTree()->CloneTree(0);

  outTree->Branch("m_svpileup", &m_svpileup, "m_svpileup/f");

  cout << "*******************************" << endl;
  cout << "*** did you fix the mistake?   " << endl;
  cout << "*** cos this file doesn't know " << endl;
  cout << "*******************************" << endl;

  for (Int_t i=0; i<eventChain->GetEntries(); i++) {
    eventChain->GetEntry(i);
    // for TTBAR
    //eventWeight=1.24*341*1000;
    // for HH
    eventWeight=2.92;
    //cout << eventWeight << ", ";
    eventWeight/=float(nevents);
    //cout << eventWeight << ", ";
    eventWeight*=0.65*0.65*btag(ptB1, etaB1)*btag(ptB2, etaB2);
    //cout << eventWeight << endl;
    //m_svpileup=mTT;
    
    outTree->Fill();

  }  
  outFile->Write();
  outFile->Save();
  
}

Double_t btag(Double_t pt, Double_t eta) {

  return ( (pt <= 20.0) * (0.000) +
           (abs(eta) <= 1.8) * (pt > 20.0 && pt <= 30) * (0.536) +
           (abs(eta) <= 1.8) * (pt > 30.0 && pt <= 40) * (0.6439) +
           (abs(eta) <= 1.8) * (pt > 40.0 && pt <= 50) * (0.6504) +
           (abs(eta) <= 1.8) * (pt > 50.0 && pt <= 60) * (0.6716) +
           (abs(eta) <= 1.8) * (pt > 60.0 && pt <= 70) * (0.6841) +
           (abs(eta) <= 1.8) * (pt > 70.0 && pt <= 80) * (0.6896) +
           (abs(eta) <= 1.8) * (pt > 80.0 && pt <= 90) * (0.6916) +
           (abs(eta) <= 1.8) * (pt > 90.0 && pt <= 100) * (0.6882) +
           (abs(eta) <= 1.8) * (pt > 100.0 && pt <= 120) * (0.6838) +
           (abs(eta) <= 1.8) * (pt > 120.0 && pt <= 140) * (0.6715) +
           (abs(eta) <= 1.8) * (pt > 140.0 && pt <= 160) * (0.6554) +
           (abs(eta) <= 1.8) * (pt > 160.0 && pt <= 180) * (0.6366) +
           (abs(eta) <= 1.8) * (pt > 180.0 && pt <= 200) * (0.6192) +
           (abs(eta) <= 1.8) * (pt > 200.0 && pt <= 250) * (0.595) +
           (abs(eta) <= 1.8) * (pt > 250.0 && pt <= 300) * (0.5551) +
           (abs(eta) <= 1.8) * (pt > 300.0 && pt <= 350) * (0.5138) +
           (abs(eta) <= 1.8) * (pt > 350.0 && pt <= 400) * (0.4884) +
           (abs(eta) <= 1.8) * (pt > 400.0 && pt <= 500) * (0.4009) +
           (abs(eta) <= 1.8) * (pt > 500.0 && pt <= 600) * (0.3459) +
           (abs(eta) <= 1.8) * (pt > 600.0 && pt <= 700) * (0.2523) +
           (abs(eta) <= 1.8) * (pt > 700.0 && pt <= 800) * (0.2404) +
           (abs(eta) <= 1.8) * (pt > 800.0 && pt <= 1000) * (0.2198) +
           (abs(eta) <= 1.8) * (pt > 1000.0 && pt <= 1400) * (0.2263) +
           (abs(eta) <= 1.8) * (pt > 1400.0 && pt <= 2000) * (0.2614) +
           (abs(eta) <= 1.8) * (pt > 2000.0) * (0.3194) +
           (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt <= 20.0) * (0.000) +
           (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 20.0 && pt <= 30) * (0.3254) +
           (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 30.0 && pt <= 40) * (0.4339) +
           (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 40.0 && pt <= 50) * (0.4499) +
           (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 50.0 && pt <= 60) * (0.4716) +
           (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 60.0 && pt <= 70) * (0.4766) +
           (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 70.0 && pt <= 80) * (0.4788) +
           (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 80.0 && pt <= 90) * (0.4863) +
           (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 90.0 && pt <= 100) * (0.4891) +
           (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 100.0 && pt <= 120) * (0.462) +
           (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 120.0 && pt <= 140) * (0.4583) +
           (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 140.0 && pt <= 160) * (0.4247) +
           (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 160.0 && pt <= 180) * (0.3775) +
           (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 180.0 && pt <= 200) * (0.3734) +
           (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 200.0 && pt <= 250) * (0.3348) +
           (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 250.0 && pt <= 300) * (0.2939) +
           (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 300.0 && pt <= 350) * (0.285) +
           (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 350.0 && pt <= 400) * (0.2421) +
           (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 400.0 && pt <= 500) * (0.1565) +
           (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 500.0 && pt <= 600) * (0.1522) +
           (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 600.0 && pt <= 700) * (0.1231) +
           (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 700.0 && pt <= 800) * (0.1607) +
           (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 800.0 && pt <= 1000) * (0.1323) +
           (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1000.0 && pt <= 1400) * (0.0) +
           (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1400.0 && pt <= 2000) * (0.0) +
           (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 2000.0) * (0.0) +
           (abs(eta) > 2.4) * (0.000));
}
