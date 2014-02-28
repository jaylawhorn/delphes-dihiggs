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

#include "modules/Delphes.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "classes/DelphesClasses.h"

#include "../MitStyleRemix.hh"
#include "../CPlot.hh"
#endif

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;

void confParse(const TString conf, vector<TString> &sampleNames, vector<TString> &sampleTitles, vector<Int_t> &sampleColors);

void newCutFlow(const TString conf="new.conf", const Int_t nSamples=4) {

  // define kinematic/plotting constants
  const Float_t TAU_PT_MIN_H = 30;
  const Float_t TAU_PT_MIN_L = 20;
  const Float_t B_PT_MIN = 30;

  const Float_t TAU_ETA_MAX = 2.5;
  const Float_t B_ETA_MAX = 2.5;
  const Int_t ETA_BINS = 12;

  const Float_t TAU_PT_MAX = 500;
  const Float_t B_PT_MAX = 500;

  enum { hadron=1, electron, muon };

  const Float_t HTT_MAX = 180;
  const Float_t HTT_MIN = 60;
  const Float_t HBB_MAX = 160;
  const Float_t HBB_MIN = 80;

  vector<TString> sampleNames;
  vector<TString> sampleTitles;
  vector<Int_t> sampleColors;

  confParse(conf, sampleNames, sampleTitles, sampleColors);

  cout << endl;
  cout << " --- Applied Cuts: --- " << endl;
  cout << "b   pT > " << B_PT_MIN << " and |eta| < " << B_ETA_MAX << endl;
  cout << "tau pT > " << TAU_PT_MIN_L << " (leptonic) or " << TAU_PT_MIN_H << " (hadronic) and |eta| < " << TAU_ETA_MAX << endl;
  cout << endl;

  UInt_t nEvents=0;

  Float_t eventWeight=1;
  Float_t met, metPhi;
  UInt_t bTag1, bTag2;
  UInt_t tauCat1, tauCat2;
  LorentzVector *recoB1=0, *recoB2=0;
  LorentzVector *recoTau1=0, *recoTau2=0;
  LorentzVector *recoLeadJet=0, *recoExtraJet=0;

  TFile *infile;
  TTree *intree;

  Double_t noCuts=0, oldCuts=0;
  Int_t iNoCuts=0, iOldCuts=0;
  Double_t jetjet=0, jetmu=0, jetele=0, elemu=0;
  Int_t iJetjet=0, iJetmu=0, iJetele=0, iElemu=0;
  Double_t extraJet=0, leadJet=0;
  Int_t iExtraJet=0, iLeadJet=0;

  Float_t bPt1=0, bPt2=0;
  Float_t tauPt1=0, tauPt2=0;

  for (UInt_t iSamp=0; iSamp<nSamples; iSamp++) { // sample loop

    TString infilename = sampleNames[iSamp];
    //cout << "Processing  " << infilename << " ..." << endl;
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

    noCuts=0; iNoCuts=0; oldCuts=0; iOldCuts=0;
    jetjet=0; jetmu=0; jetele=0; elemu=0;
    iJetjet=0; iJetmu=0; iJetele=0; iElemu=0;
    extraJet=0; iExtraJet=0; leadJet=0; iLeadJet=0;

    for(UInt_t iEntry=0; iEntry<intree->GetEntries(); iEntry++) { // entry loop
      intree->GetEntry(iEntry);

      noCuts+=eventWeight;
      iNoCuts++;

      // skip events that don't have 2 reco b's and 2 reco tau's
      if ( (recoB1->Pt()==999) || (recoB2->Pt()==999) ) continue;
      if ( (recoTau1->Pt()==999) || (recoTau2->Pt()==999) ) continue;

      bPt1=recoB1->Pt(); 
      bPt2=recoB2->Pt();
      tauPt1 = recoTau1->Pt();
      tauPt2 = recoTau2->Pt();

      if ( ( fabs( recoTau1->Eta() ) > TAU_ETA_MAX ) || ( fabs(recoTau2->Eta() ) > TAU_ETA_MAX ) ) continue;
      if ( ( fabs( recoB1->Eta() ) > B_ETA_MAX ) || ( fabs( recoB2->Eta() ) > B_ETA_MAX ) ) continue;

      if ( ( bPt1 < B_PT_MIN ) || ( bPt2 < B_PT_MIN ) ) continue;

      //if ( (tauPt1 > TAU_PT_MIN_H) && (tauPt2 > TAU_PT_MIN_H) ) { oldCuts+=eventWeight; iOldCuts++; }

      // jet-jet
      if ( (tauCat1==hadron) && (tauCat2==hadron) ) {
	if ( (tauPt1 > TAU_PT_MIN_H) && (tauPt2 > TAU_PT_MIN_H) ) { iJetjet++; jetjet+=eventWeight; }
      }
      // jet-mu
      else if ( (tauCat1==hadron) && (tauCat2==muon) ) {
	if ( (tauPt1 > TAU_PT_MIN_H) && (tauPt2 > TAU_PT_MIN_L) ) { iJetmu++; jetmu+=eventWeight; }
      }
      else if ( (tauCat1==muon) && (tauCat2==hadron) ) {
	if ( (tauPt1 > TAU_PT_MIN_L) && (tauPt2 > TAU_PT_MIN_H) ) { iJetmu++; jetmu+=eventWeight; }
      }
      // jet-ele
      else if ( (tauCat1==hadron) && (tauCat2==electron) ) {
	if ( (tauPt1 > TAU_PT_MIN_H) && (tauPt2 > TAU_PT_MIN_L) ) { iJetele++; jetele+=eventWeight; }
      }
      else if ( (tauCat1==electron) && (tauCat2==hadron) ) {
	if ( (tauPt1 > TAU_PT_MIN_L) && (tauPt2 > TAU_PT_MIN_H) ) { iJetele++; jetele+=eventWeight; }
      }
      // ele-mu
      else if ( (tauCat1==muon) && (tauCat2==electron) ) {
	if ( (tauPt1 > TAU_PT_MIN_L) && (tauPt2 > TAU_PT_MIN_L) ) { iElemu++; elemu+=eventWeight; }
      }
      else if ( (tauCat1==electron) && (tauCat2==muon) ) {
	if ( (tauPt1 > TAU_PT_MIN_L) && (tauPt2 > TAU_PT_MIN_L) ) { iElemu++; elemu+=eventWeight; }
      }
      else continue;
      //if ( (tauCat1 == muon) && (tauCat2 == muon) ) continue;
      //if ( (tauCat1 == electron) && (tauCat2 == electron) ) continue;

      if ((recoExtraJet->Pt()!=999) && (recoExtraJet->Pt()>30)) { extraJet+=eventWeight; iExtraJet++; }

      if ((recoExtraJet->Pt()!=999) && (recoExtraJet->Pt()>50)) { leadJet+=eventWeight; iLeadJet++; }

    } // end entry loop

    cout << "-------" << endl;
    cout << sampleTitles[iSamp] << " event yields at 3000/fb: "<< endl;
    //cout << "total (no cuts):         " << noCuts*3000 << endl;
    //cout << "old cuts:                " << oldCuts*3000 << endl;
    cout << "Tau-tagged Dijet:        " << jetjet*3000 << endl;
    cout << "Tau-tagged jet-muon:     " << jetmu*3000 << endl;
    cout << "Tau-tagged jet-electron: " << jetele*3000 << endl;
    cout << "Electron-muon:           " << elemu*3000 << endl;
    cout << "All four channels:       " << 3000*(jetjet+jetmu+jetele+elemu) << endl;
    cout << "+an extra jet (pT>30):   " << 3000*extraJet << endl;
    cout << "              (pT>50):   " << 3000*leadJet << endl;

    delete infile;
    infile=0, intree=0;

  } // end sample loop

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
