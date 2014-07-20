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

//#include "../Utils/ttMVA.h"

#endif

void combinefiles(const TString input="temp.txt",
		  const TString outputfile="test.root") {

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
  Int_t nCentral=0, nBtag=0; Int_t centB=0;

  Int_t tauCat1=0, tauCat2=0;
  Int_t bTag1=0, bTag2=0, bTag3=0, bTag4=0;
  Int_t jbTag_tt1=0, jbTag_tt2=0, jbTag_6j1=0, jbTag_6j2=0;

  Double_t mt2=0; Double_t mt2puppi=0; Double_t mt2pileup=0;
  Double_t m_sv=0; Double_t m_svpileup=0; Double_t m_svpuppi=0;
  
  Float_t ptTau1, ptTau2, ptB1, ptB2, ptB3, ptB4, ptG1, ptG2, ptJet_tt1, ptJet_tt2, ptJet_6j1, ptJet_6j2, tauIso1, tauIso2;
  Float_t etaTau1, etaTau2, etaB1, etaB2, etaB3, etaB4, etaG1, etaG2, etaJet_tt1, etaJet_tt2, etaJet_6j1, etaJet_6j2;
  Float_t phiTau1, phiTau2, phiB1, phiB2, phiB3, phiB4, phiG1, phiG2, phiJet_tt1, phiJet_tt2, phiJet_6j1, phiJet_6j2;
  Float_t mTau1, mTau2, mB1, mB2, mB3, mB4, eG1, eG2, mJet_tt1, mJet_tt2, mJet_6j1, mJet_6j2;

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
  Float_t bdtVal=0;

  TFile *outFile = new TFile(outputfile, "RECREATE");

  countChain->SetBranchAddress("n", &n);

  Int_t nevents=0;

  for (Int_t i=0; i<countChain->GetEntries(); i++) {
    countChain->GetEntry(i);
    nevents+=n;
  }
  cout << nevents << endl;

  eventChain->SetBranchAddress("eventWeight",    &eventWeight);
  eventChain->SetBranchAddress("sampleNo",       &sampleNo);
  eventChain->SetBranchAddress("genInfo",        &genInfo);
  eventChain->SetBranchAddress("isBBTT",         &isBBTT);
  eventChain->SetBranchAddress("isBBGG",         &isBBGG);
  eventChain->SetBranchAddress("isBBBB",         &isBBBB);
  eventChain->SetBranchAddress("isVBFTT",        &isVBFTT);
  eventChain->SetBranchAddress("isVBF4B",        &isVBF4B);
  eventChain->SetBranchAddress("eventType",      &eventType);

  eventChain->SetBranchAddress("met",            &met);
  eventChain->SetBranchAddress("metPhi",         &metPhi);
  eventChain->SetBranchAddress("pileupmet",      &pileupMet);
  eventChain->SetBranchAddress("pileupmetPhi",   &pileupMetPhi);
  eventChain->SetBranchAddress("puppiMet",       &ppMet);
  eventChain->SetBranchAddress("puppiMetPhi",    &ppMetPhi);

  eventChain->SetBranchAddress("ptTau1",         &ptTau1);
  eventChain->SetBranchAddress("etaTau1",        &etaTau1);
  eventChain->SetBranchAddress("phiTau1",        &phiTau1);
  eventChain->SetBranchAddress("mTau1",          &mTau1);
  eventChain->SetBranchAddress("tauCat1",        &tauCat1);
  eventChain->SetBranchAddress("tauIso1",        &tauIso1);
  eventChain->SetBranchAddress("ptTau2",         &ptTau2);
  eventChain->SetBranchAddress("etaTau2",        &etaTau2);
  eventChain->SetBranchAddress("phiTau2",        &phiTau2);
  eventChain->SetBranchAddress("mTau2",          &mTau2);
  eventChain->SetBranchAddress("tauCat2",        &tauCat2);
  eventChain->SetBranchAddress("tauIso2",        &tauIso2);

  eventChain->SetBranchAddress("ptG1",           &ptG1);
  eventChain->SetBranchAddress("etaG1",          &etaG1);
  eventChain->SetBranchAddress("phiG1",          &phiG1);
  eventChain->SetBranchAddress("eG1",            &eG1);
  eventChain->SetBranchAddress("gamIso1",        &gamIso1);

  eventChain->SetBranchAddress("ptG2",           &ptG2);
  eventChain->SetBranchAddress("etaG2",          &etaG2);
  eventChain->SetBranchAddress("phiG2",          &phiG2);
  eventChain->SetBranchAddress("eG2",            &eG2);
  eventChain->SetBranchAddress("gamIso2",        &gamIso2);

  eventChain->SetBranchAddress("ptB1",           &ptB1);
  eventChain->SetBranchAddress("etaB1",          &etaB1);
  eventChain->SetBranchAddress("phiB1",          &phiB1);
  eventChain->SetBranchAddress("mB1",            &mB1);
  eventChain->SetBranchAddress("bTag1",          &bTag1);

  eventChain->SetBranchAddress("ptB2",           &ptB2);
  eventChain->SetBranchAddress("etaB2",          &etaB2);
  eventChain->SetBranchAddress("phiB2",          &phiB2);
  eventChain->SetBranchAddress("mB2",            &mB2);
  eventChain->SetBranchAddress("bTag2",          &bTag2);

  eventChain->SetBranchAddress("ptB3",           &ptB3);
  eventChain->SetBranchAddress("etaB3",          &etaB3);
  eventChain->SetBranchAddress("phiB3",          &phiB3);
  eventChain->SetBranchAddress("mB3",            &mB3);
  eventChain->SetBranchAddress("bTag3",          &bTag3);

  eventChain->SetBranchAddress("ptB4",           &ptB4);
  eventChain->SetBranchAddress("etaB4",          &etaB4);
  eventChain->SetBranchAddress("phiB4",          &phiB4);
  eventChain->SetBranchAddress("mB4",            &mB4);
  eventChain->SetBranchAddress("bTag4",          &bTag4);

  eventChain->SetBranchAddress("ptJet_tt1",      &ptJet_tt1);
  eventChain->SetBranchAddress("etaJet_tt1",     &etaJet_tt1);
  eventChain->SetBranchAddress("phiJet_tt1",     &phiJet_tt1);
  eventChain->SetBranchAddress("mJet_tt1",       &mJet_tt1);
  eventChain->SetBranchAddress("jbTag_tt1",      &jbTag_tt1);

  eventChain->SetBranchAddress("ptJet_tt2",      &ptJet_tt2);
  eventChain->SetBranchAddress("etaJet_tt2",     &etaJet_tt2);
  eventChain->SetBranchAddress("phiJet_tt2",     &phiJet_tt2);
  eventChain->SetBranchAddress("mJet_tt2",       &mJet_tt2);
  eventChain->SetBranchAddress("jbTag_tt2",      &jbTag_tt2);

  eventChain->SetBranchAddress("ptJet_6j1",      &ptJet_6j1);
  eventChain->SetBranchAddress("etaJet_6j1",     &etaJet_6j1);
  eventChain->SetBranchAddress("phiJet_6j1",     &phiJet_6j1);
  eventChain->SetBranchAddress("mJet_6j1",       &mJet_6j1);
  eventChain->SetBranchAddress("jbTag_6j1",      &jbTag_6j1);

  eventChain->SetBranchAddress("ptJet_6j2",      &ptJet_6j2);
  eventChain->SetBranchAddress("etaJet_6j2",     &etaJet_6j2);
  eventChain->SetBranchAddress("phiJet_6j2",     &phiJet_6j2);
  eventChain->SetBranchAddress("mJet_6j2",       &mJet_6j2);
  eventChain->SetBranchAddress("jbTag_6j2",      &jbTag_6j2);

  eventChain->SetBranchAddress("ptTau1_gen",     &ptTau1_gen);
  eventChain->SetBranchAddress("etaTau1_gen",    &etaTau1_gen);
  eventChain->SetBranchAddress("phiTau1_gen",    &phiTau1_gen);
  eventChain->SetBranchAddress("mTau1_gen",      &mTau1_gen);

  eventChain->SetBranchAddress("ptTau2_gen",     &ptTau2_gen);
  eventChain->SetBranchAddress("etaTau2_gen",    &etaTau2_gen);
  eventChain->SetBranchAddress("phiTau2_gen",    &phiTau2_gen);
  eventChain->SetBranchAddress("mTau2_gen",      &mTau2_gen);

  eventChain->SetBranchAddress("ptTau1_genJet",  &ptTau1_genJet);
  eventChain->SetBranchAddress("etaTau1_genJet", &etaTau1_genJet);
  eventChain->SetBranchAddress("phiTau1_genJet", &phiTau1_genJet);
  eventChain->SetBranchAddress("mTau1_genJet",   &mTau1_genJet);

  eventChain->SetBranchAddress("ptTau2_genJet",  &ptTau2_genJet);
  eventChain->SetBranchAddress("etaTau2_genJet", &etaTau2_genJet);
  eventChain->SetBranchAddress("phiTau2_genJet", &phiTau2_genJet);
  eventChain->SetBranchAddress("mTau2_genJet",   &mTau2_genJet);

  eventChain->SetBranchAddress("ptG1_gen",       &ptG1_gen);
  eventChain->SetBranchAddress("etaG1_gen",      &etaG1_gen);
  eventChain->SetBranchAddress("phiG1_gen",      &phiG1_gen);
  eventChain->SetBranchAddress("eG1_gen",        &eG1_gen);

  eventChain->SetBranchAddress("ptG2_gen",       &ptG2_gen);
  eventChain->SetBranchAddress("etaG2_gen",      &etaG2_gen);
  eventChain->SetBranchAddress("phiG2_gen",      &phiG2_gen);
  eventChain->SetBranchAddress("eG2_gen",        &eG2_gen);

  eventChain->SetBranchAddress("ptB1_gen",       &ptB1_gen);
  eventChain->SetBranchAddress("etaB1_gen",      &etaB1_gen);
  eventChain->SetBranchAddress("phiB1_gen",      &phiB1_gen);
  eventChain->SetBranchAddress("mB1_gen",        &mB1_gen);
  eventChain->SetBranchAddress("iHmatch1",       &iHmatch1);

  eventChain->SetBranchAddress("ptB2_gen",       &ptB2_gen);
  eventChain->SetBranchAddress("etaB2_gen",      &etaB2_gen);
  eventChain->SetBranchAddress("phiB2_gen",      &phiB2_gen);
  eventChain->SetBranchAddress("mB2_gen",        &mB2_gen);
  eventChain->SetBranchAddress("iHmatch2",       &iHmatch2);

  eventChain->SetBranchAddress("ptB3_gen",       &ptB3_gen);
  eventChain->SetBranchAddress("etaB3_gen",      &etaB3_gen);
  eventChain->SetBranchAddress("phiB3_gen",      &phiB3_gen);
  eventChain->SetBranchAddress("mB3_gen",        &mB3_gen);
  eventChain->SetBranchAddress("iHmatch3",       &iHmatch3);

  eventChain->SetBranchAddress("ptB4_gen",       &ptB4_gen);
  eventChain->SetBranchAddress("etaB4_gen",      &etaB4_gen);
  eventChain->SetBranchAddress("phiB4_gen",      &phiB4_gen);
  eventChain->SetBranchAddress("mB4_gen",        &mB4_gen);
  eventChain->SetBranchAddress("iHmatch4",       &iHmatch4);

  eventChain->SetBranchAddress("ptH1_gen",       &ptH1_gen);
  eventChain->SetBranchAddress("etaH1_gen",      &etaH1_gen);
  eventChain->SetBranchAddress("phiH1_gen",      &phiH1_gen);
  eventChain->SetBranchAddress("mH1_gen",        &mH1_gen);

  eventChain->SetBranchAddress("ptH2_gen",       &ptH2_gen);
  eventChain->SetBranchAddress("etaH2_gen",      &etaH2_gen);
  eventChain->SetBranchAddress("phiH2_gen",      &phiH2_gen);
  eventChain->SetBranchAddress("mH2_gen",        &mH2_gen);

  eventChain->SetBranchAddress("ptJet_tt1_gen",  &ptJet_tt1_gen);
  eventChain->SetBranchAddress("etaJet_tt1_gen", &etaJet_tt1_gen);
  eventChain->SetBranchAddress("phiJet_tt1_gen", &phiJet_tt1_gen);
  eventChain->SetBranchAddress("mJet_tt1_gen",   &mJet_tt1_gen);

  eventChain->SetBranchAddress("ptJet_tt2_gen",  &ptJet_tt2_gen);
  eventChain->SetBranchAddress("etaJet_tt2_gen", &etaJet_tt2_gen);
  eventChain->SetBranchAddress("phiJet_tt2_gen", &phiJet_tt2_gen);
  eventChain->SetBranchAddress("mJet_tt2_gen",   &mJet_tt2_gen);

  eventChain->SetBranchAddress("ptJet_6j1_gen",  &ptJet_6j1_gen);
  eventChain->SetBranchAddress("etaJet_6j1_gen", &etaJet_6j1_gen);
  eventChain->SetBranchAddress("phiJet_6j1_gen", &phiJet_6j1_gen);
  eventChain->SetBranchAddress("mJet_6j1_gen",   &mJet_6j1_gen);

  eventChain->SetBranchAddress("ptJet_6j2_gen",  &ptJet_6j2_gen);
  eventChain->SetBranchAddress("etaJet_6j2_gen", &etaJet_6j2_gen);
  eventChain->SetBranchAddress("phiJet_6j2_gen", &phiJet_6j2_gen);
  eventChain->SetBranchAddress("mJet_6j2_gen",   &mJet_6j2_gen);

  eventChain->SetBranchAddress("ptTT",           &ptTT);
  eventChain->SetBranchAddress("etaTT",          &etaTT);
  eventChain->SetBranchAddress("phiTT",          &phiTT);
  eventChain->SetBranchAddress("mTT",            &mTT);

  eventChain->SetBranchAddress("ptBB1",          &ptBB1);
  eventChain->SetBranchAddress("etaBB1",         &etaBB1);
  eventChain->SetBranchAddress("phiBB1",         &phiBB1);
  eventChain->SetBranchAddress("mBB1",           &mBB1);

  eventChain->SetBranchAddress("ptBB2",          &ptBB2);
  eventChain->SetBranchAddress("etaBB2",         &etaBB2);
  eventChain->SetBranchAddress("phiBB2",         &phiBB2);
  eventChain->SetBranchAddress("mBB2",           &mBB2);

  eventChain->SetBranchAddress("ptGG",           &ptGG);
  eventChain->SetBranchAddress("etaGG",          &etaGG);
  eventChain->SetBranchAddress("phiGG",          &phiGG);
  eventChain->SetBranchAddress("mGG",            &mGG);

  eventChain->SetBranchAddress("ptJJ_tt",        &ptJJ_tt);
  eventChain->SetBranchAddress("etaJJ_tt",       &etaJJ_tt);
  eventChain->SetBranchAddress("phiJJ_tt",       &phiJJ_tt);
  eventChain->SetBranchAddress("mJJ_tt",         &mJJ_tt);

  eventChain->SetBranchAddress("ptJJ_6j",        &ptJJ_6j);
  eventChain->SetBranchAddress("etaJJ_6j",       &etaJJ_6j);
  eventChain->SetBranchAddress("phiJJ_6j",       &phiJJ_6j);
  eventChain->SetBranchAddress("mJJ_6j",         &mJJ_6j);

  eventChain->SetBranchAddress("ptHH",           &ptHH);
  eventChain->SetBranchAddress("etaHH",          &etaHH);
  eventChain->SetBranchAddress("phiHH",          &phiHH);
  eventChain->SetBranchAddress("mHH",            &mHH);

  eventChain->SetBranchAddress("mt2",            &mt2);
  eventChain->SetBranchAddress("mt2puppi",       &mt2puppi);
  eventChain->SetBranchAddress("mt2pileup",      &mt2pileup);
  
  eventChain->SetBranchAddress("m_sv",           &m_sv);
  eventChain->SetBranchAddress("m_svpileup",     &m_svpileup);
  eventChain->SetBranchAddress("m_svpuppi",      &m_svpuppi);

  eventChain->SetBranchAddress("nBtag",          &nBtag);
  eventChain->SetBranchAddress("nCentral",       &nCentral);
  eventChain->SetBranchAddress("centB",          &centB);
  eventChain->SetBranchAddress("dEta_tt",        &dEta_tt);
  eventChain->SetBranchAddress("dEta_6j",        &dEta_6j);
  
  eventChain->GetEntry(0);

  TTree *outTree=(TTree*)eventChain->GetTree()->CloneTree(0);

  outTree->Branch("bdtVal", &bdtVal, "bdtVal/f");

  /*  ttMVA::MVAType why=ttMVA::kMuTau;

  ttMVA *testMVA = new ttMVA();
  testMVA->Intialize(why);*/

  for (Int_t i=0; i<eventChain->GetEntries(); i++) {
    eventChain->GetEntry(i);

    if (isBBTT==0) continue;

    eventWeight/=float(nevents);
    bdtVal=999;
    //bdtVal=testMVA->GetBDTValue(mTT, mBB1, mHH, mt2, met, m_sv, ptHH, TMath::Sqrt((etaB1-etaB2)*(etaB1-etaB2)+(phiB1-phiB2)*(phiB1-phiB2)), TMath::Sqrt((etaB1-etaB2)*(etaB1-etaB2)+(phiB1-phiB2)*(phiB1-phiB2)));
    //cout << bdtVal << endl;

    outTree->Fill();
  }

  outFile->Write();
  outFile->Save();

}
