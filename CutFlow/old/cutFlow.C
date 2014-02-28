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

void cutFlow(const TString conf="new.conf", const Int_t nSamples=4) {

  // define kinematic/plotting constants
  const Float_t TAU_PT_MIN_H = 30;
  const Float_t TAU_PT_MIN_L = 30;
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
  UInt_t bTag1, bTag2;
  UInt_t tauDecayCat1, tauDecayCat2;
  LorentzVector *genB1=0, *genB2=0, *recoB1=0, *recoB2=0;
  LorentzVector *genTau1=0, *genTau2=0, *genDecayTau1=0, *genDecayTau2=0, *recoTau1=0, *recoTau2=0;
  LorentzVector *tauHiggs=0, *bHiggs=0;

  TFile *infile;
  TTree *intree;

  Double_t hadDecay=0, eleDecay=0, muonDecay=0;
  Double_t noCuts=0, twoProd=0, fourProd=0;
  Double_t bEta=0, bPt=0, tauEta=0, tauPt=0;
  Double_t ttMass=0, bbMass=0; 

  Int_t inoCuts=0, itwoProd=0, ifourProd=0;
  Int_t ibEta=0, ibPt=0, itauEta=0, itauPt=0;
  Int_t ittMass=0, ibbMass=0; 

  Float_t bCorPt1=0, bCorPt2=0;
  Float_t tauCorPt1=0, tauCorPt2=0;

  //TCanvas *c = MakeCanvas("c", "c", 1200, 600);

  TH1D *hTTMass = new TH1D("hTTMass", "hTTMass", 200, 0, 400);
  TH1D *hBBMass = new TH1D("hBBMass", "hBBMass", 200, 0, 400);

  for (UInt_t iSamp=0; iSamp<nSamples; iSamp++) { // sample loop

    TString infilename = sampleNames[iSamp];
    //cout << "Processing  " << infilename << " ..." << endl;
    infile = new TFile(infilename); assert(infile);
    intree = (TTree*) infile->Get("Events"); assert(intree);
 
    intree->SetBranchAddress("eventWeight",    &eventWeight);
    intree->SetBranchAddress("bTag1",          &bTag1);
    intree->SetBranchAddress("bTag2",          &bTag2);
    intree->SetBranchAddress("genB1",          &genB1);
    intree->SetBranchAddress("genB2",          &genB2);
    intree->SetBranchAddress("recoB1",         &recoB1);
    intree->SetBranchAddress("recoB2",         &recoB2);
    intree->SetBranchAddress("tauDecayCat1",   &tauDecayCat1);
    intree->SetBranchAddress("tauDecayCat2",   &tauDecayCat2);
    intree->SetBranchAddress("genTau1",        &genTau1);
    intree->SetBranchAddress("genTau2",        &genTau2);
    intree->SetBranchAddress("genDecayTau1",   &genDecayTau1);
    intree->SetBranchAddress("genDecayTau2",   &genDecayTau2);
    intree->SetBranchAddress("recoTau1",       &recoTau1);
    intree->SetBranchAddress("recoTau2",       &recoTau2);

    noCuts=0; twoProd=0; fourProd=0;
    bEta=0; bPt=0; tauEta=0; tauPt=0;
    ttMass=0; bbMass=0;

    inoCuts=0; itwoProd=0; ifourProd=0;
    ibEta=0; ibPt=0; itauEta=0; itauPt=0;
    ittMass=0; ibbMass=0;

    for(UInt_t iEntry=0; iEntry<intree->GetEntries(); iEntry++) { // entry loop
      intree->GetEntry(iEntry);

      noCuts+=eventWeight;
      inoCuts++;

      // skip events that don't have 2 reco b's and 2 reco tau's
      if ( (recoB1->Pt()==999) || (recoB2->Pt()==999) ) continue;
      if ( (tauDecayCat1 == muon) && (tauDecayCat2 == muon) ) continue;
      if ( (tauDecayCat1 == electron) && (tauDecayCat2 == electron) ) continue;

      twoProd+=eventWeight;
      itwoProd++;

      if ( (recoTau1->Pt()==999) || (recoTau2->Pt()==999) ) continue;

      fourProd+=eventWeight;
      ifourProd++;

      bCorPt1=recoB1->Pt(); 
      bCorPt2=recoB2->Pt();
      tauCorPt1 = recoTau1->Pt();
      tauCorPt2 = recoTau2->Pt();

      if ( ( fabs( recoTau1->Eta() ) > TAU_ETA_MAX ) || ( fabs(recoTau2->Eta() ) > TAU_ETA_MAX ) ) continue;

      tauEta+=eventWeight;
      itauEta++;

      if (( tauDecayCat1 == hadron ) && (recoTau1->Pt() < TAU_PT_MIN_H )) continue;
      else if (( tauDecayCat2 == hadron ) && (recoTau2->Pt() < TAU_PT_MIN_H )) continue;
      else if ( ( recoTau1->Pt() < TAU_PT_MIN_L ) || ( recoTau2->Pt() < TAU_PT_MIN_L ) ) continue;

      tauPt+=eventWeight;
      itauPt++;

      if ( ( fabs( recoB1->Eta() ) > B_ETA_MAX ) || ( fabs( recoB2->Eta() ) > B_ETA_MAX ) ) continue;

      bEta+=eventWeight;
      ibEta++;

      if ( ( bCorPt1 < B_PT_MIN ) || ( bCorPt2 < B_PT_MIN ) ) continue;

      bPt+=eventWeight;
      ibPt++;

      LorentzVector vTau1(tauCorPt1, recoTau1->Eta(), recoTau1->Phi(), recoTau1->M());
      LorentzVector vTau2(tauCorPt2, recoTau2->Eta(), recoTau2->Phi(), recoTau2->M());
      LorentzVector vTauHiggs = vTau1 + vTau2;

      hTTMass->Fill(vTauHiggs.M());

      LorentzVector vB1(bCorPt1, recoB1->Eta(), recoB1->Phi(), recoB1->M());
      LorentzVector vB2(bCorPt2, recoB2->Eta(), recoB2->Phi(), recoB2->M());
      LorentzVector vBHiggs = vB1 + vB2;

      hBBMass->Fill(vBHiggs.M());

      if ( ( vTauHiggs.M() > HTT_MAX ) || ( vTauHiggs.M() < HTT_MIN ) ) continue;

      ttMass+=eventWeight;
      ittMass++;

      if ( ( vBHiggs.M() > HBB_MAX ) || ( vBHiggs.M() < HBB_MIN ) ) continue;

      bbMass+=eventWeight;
      ibbMass++;      

    } // end entry loop
    /*
    c->Divide(2,1);

    c->cd(1)->SetPad(0,0,0.5,1);
    c->cd(1)->SetTopMargin(0.1);
    c->cd(1)->SetBottomMargin(0.2);
    c->cd(1)->SetLeftMargin(0.2);
    c->cd(1)->SetRightMargin(0.1);
    c->cd(2)->SetPad(0.5,0,1,1);
    c->cd(2)->SetTopMargin(0.1);
    c->cd(2)->SetBottomMargin(0.2);
    c->cd(2)->SetLeftMargin(0.15);
    c->cd(2)->SetRightMargin(0.15);

    c->cd(1);

    TLine *ttMin = new TLine(HTT_MIN, 0, HTT_MIN, 25);
    ttMin->SetLineColor(kRed);
    ttMin->SetLineWidth(2);
    TLine *ttMax = new TLine(HTT_MAX, 0, HTT_MAX, 25);
    ttMax->SetLineColor(kRed);
    ttMax->SetLineWidth(2);

    hTTMass->SetTitle("");
    hTTMass->GetXaxis()->SetTitle("M(#tau#tau) [GeV/c^{2}]");
    hTTMass->GetYaxis()->SetTitle("Events");

    hTTMass->Draw();
    ttMin->Draw("same");
    ttMax->Draw("same");

    c->cd(2);

    TLine *bbMin = new TLine(HBB_MIN, 0, HBB_MIN, 31);
    bbMin->SetLineColor(kRed);
    bbMin->SetLineWidth(2);
    TLine *bbMax = new TLine(HBB_MAX, 0, HBB_MAX, 31);
    bbMax->SetLineColor(kRed);
    bbMax->SetLineWidth(2);

    hBBMass->SetTitle("");
    hBBMass->GetXaxis()->SetTitle("M(bb) [GeV/c^{2}]");

    hBBMass->Draw();
    bbMin->Draw("same");
    bbMax->Draw("same");

    c->SaveAs("mass_dists.png");
    */    
    cout << sampleTitles[iSamp].Data() <<  endl;
    //cout << "                          " << "Scaled \tUnscaled" << endl; 
    cout << "Events w/ bb, tautau    : " << 3000*noCuts << "\t\t" << inoCuts << endl;
    //cout << " ... after tau eta cut  : " << tauEta << "\t" << itauEta << endl;
    cout << " ... after tau pt cut   : " << 3000*tauPt << "\t" << itauPt << endl;
    //cout << " ... after b eta cut    : " << bEta << "\t" << ibEta << endl;
    cout << " ... after b pt cut     : " << 3000*bPt << "\t" << ibPt << endl;
    //cout << " ... after tt mass cut  : " << ttMass << "\t" << ittMass << endl;
    //cout << " ... after bb mass cut  : " << bbMass << "\t" << ibbMass << endl;

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
