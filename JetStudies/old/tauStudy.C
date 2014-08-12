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
#include "Math/LorentzVector.h"

#include "modules/Delphes.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "classes/DelphesClasses.h"

#include "MitStyleRemix.hh"
#include "bJetScaleCorr.h"
#endif

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;

void confParse(const TString conf, vector<TString> &sampleNames, vector<TString> &sampleTitles, vector<Int_t> &sampleColors);

void tauStudy(const TString conf="new.conf") {

  // tau decay modes
  enum { hadron=1, electron, muon };

  const Int_t nSamples=1;
  
  vector<TString> sampleNames;
  vector<TString> sampleTitles;
  vector<Int_t> sampleColors;

  confParse(conf, sampleNames, sampleTitles, sampleColors);

  TProfile *hEffPt[nSamples];
  TProfile *hEffEta[nSamples];

  TProfile *hJetResPt[nSamples];
  TProfile *hJetResEta[nSamples];

  TH1D *hPt[nSamples];

  char hname[100];

  // define kinematic/plotting constants
  const Float_t PT_MAX  = 300;
  const Float_t PT_MIN  = 0;
  const Int_t   PT_BIN  = 150;
  const Float_t ETA_MAX = 2.5;
  const Float_t ETA_MIN = -2.5;
  const Int_t   ETA_BIN = 16;

  Double_t jetCorr1Pt, jetCorr2Pt;

  for(UInt_t iSamp=0; iSamp<nSamples; iSamp++) {

    sprintf(hname, "hEffPt_%s", sampleTitles[iSamp].Data()); hEffPt[iSamp]= new TProfile(hname, hname, PT_BIN, PT_MIN, PT_MAX);
    sprintf(hname, "hEffEta_%s", sampleTitles[iSamp].Data()); hEffEta[iSamp]= new TProfile(hname, hname, ETA_BIN, ETA_MIN, ETA_MAX);

    sprintf(hname, "hJetResPt_%s", sampleTitles[iSamp].Data()); hJetResPt[iSamp] = new TProfile(hname, hname, PT_BIN, PT_MIN, PT_MAX);
    sprintf(hname, "hJetResEta_%s", sampleTitles[iSamp].Data()); hJetResEta[iSamp] = new TProfile(hname, hname, ETA_BIN, ETA_MIN, ETA_MAX);

    sprintf(hname, "hPt_%s", sampleTitles[iSamp].Data()); hPt[iSamp] = new TH1D(hname, hname, PT_BIN, PT_MIN, PT_MAX);

  }

  UInt_t eventNum;
  UInt_t bTag1, bTag2;
  UInt_t tauCat1, tauCat2;
  LorentzVector *genB1=0, *genB2=0, *recoB1=0, *recoB2=0;
  LorentzVector *genTau1=0, *genTau2=0, *genDecayTau1=0, *genDecayTau2=0, *recoTau1=0, *recoTau2=0;

  TFile *infile;
  TTree *intree;

  for (UInt_t iSamp=0; iSamp<nSamples; iSamp++) { // sample loop

    TString infilename = sampleNames[iSamp];
    cout << "Processing  " << infilename << " ..." << endl;
    infile = new TFile(infilename); assert(infile);
    intree = (TTree*) infile->Get("Events"); assert(intree);
 
    //intree->SetBranchAddress("eventNum",       &eventNum);
    intree->SetBranchAddress("bTag1",          &bTag1);
    intree->SetBranchAddress("bTag2",          &bTag2);
    intree->SetBranchAddress("genB1",          &genB1);
    intree->SetBranchAddress("genB2",          &genB2);
    intree->SetBranchAddress("recoB1",         &recoB1);
    intree->SetBranchAddress("recoB2",         &recoB2);
    intree->SetBranchAddress("tauCat1",        &tauCat1);
    intree->SetBranchAddress("tauCat2",        &tauCat2);
    intree->SetBranchAddress("genTau1",        &genTau1);
    intree->SetBranchAddress("genTau2",        &genTau2);
    intree->SetBranchAddress("genDecayTau1",   &genDecayTau1);
    intree->SetBranchAddress("genDecayTau2",   &genDecayTau2);
    intree->SetBranchAddress("recoTau1",       &recoTau1);
    intree->SetBranchAddress("recoTau2",       &recoTau2);

    for(UInt_t iEntry=0; iEntry<intree->GetEntries(); iEntry++) { // entry loop
      intree->GetEntry(iEntry);

      // jet resolution
      //cout << tauCat1 << " " << tauCat2 << endl;
      if ((recoTau1->Pt()!=999)) {

	//jetCorr1Pt = recoTau1->Pt()*getJetScaleFactor(recoTau1->Pt(), recoTau1->Eta());
	jetCorr1Pt=recoTau1->Pt();

	hPt[iSamp]->Fill(recoTau1->Pt());

	hJetResPt[iSamp]->Fill(genDecayTau1->Pt(),(recoTau1->Pt()-genDecayTau1->Pt())/genDecayTau1->Pt());
        hJetResEta[iSamp]->Fill(genDecayTau1->Eta(),(recoTau1->Pt()-genDecayTau1->Pt())/genDecayTau1->Pt());

      }

      if ((recoTau2->Pt()!=999)) {

	//jetCorr2Pt = recoTau2->Pt()*getJetScaleFactor(recoTau2->Pt(), recoTau2->Eta());
	jetCorr2Pt=recoTau2->Pt();

	hPt[iSamp]->Fill(recoTau2->Pt());

	hJetResPt[iSamp]->Fill(genDecayTau2->Pt(),(recoTau2->Pt()-genDecayTau2->Pt())/genDecayTau2->Pt());
        hJetResEta[iSamp]->Fill(genDecayTau2->Eta(),(recoTau2->Pt()-genDecayTau2->Pt())/genDecayTau2->Pt());

      }

      //tagging efficiency
      /*
      if ((tauDecayCat1==1) && (recoTau1->Pt()!=999)) {
	hEffPt[iSamp]->Fill(genDecayTau1->Pt(),1);
	hEffEta[iSamp]->Fill(genDecayTau1->Eta(),1);
      }
      else {
	hEffPt[iSamp]->Fill(genDecayTau1->Pt(),0);
	hEffEta[iSamp]->Fill(genDecayTau1->Eta(),0);
      }

      if ((tauDecayCat2==1) && (recoTau2->Pt()!=999)) {
	hEffPt[iSamp]->Fill(genDecayTau2->Pt(),1);
	hEffEta[iSamp]->Fill(genDecayTau2->Eta(),1);
      }
      else {
	hEffPt[iSamp]->Fill(genDecayTau2->Pt(),0);
	hEffEta[iSamp]->Fill(genDecayTau2->Eta(),0);
      }
      */
    } // end entry loop
    delete infile;
    infile=0, intree=0;

  } // end sample loop

  char pname[100];
  char xlabel[100];
  char ylabel[100];
  
  TCanvas *c = MakeCanvas("c", "c", 800, 600);
  TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
  leg->SetShadowColor(0);
  leg->SetFillColor(0);

  // tau jet resolution as a function of pt (uncorrected)
  sprintf(xlabel, "generator level tau jet P_{T}");
  sprintf(ylabel, "(reco P_{T} - gen P_{T})/gen P_{T}");
  sprintf(pname, "tauJetResPt");

  for (UInt_t iSamp=0; iSamp<nSamples; iSamp++) {

    hJetResPt[iSamp]->SetLineColor(sampleColors[iSamp]);
    hJetResPt[iSamp]->SetMarkerColor(sampleColors[iSamp]);
    hJetResPt[iSamp]->SetMarkerSize(1);
    leg->AddEntry(hJetResPt[iSamp], sampleTitles[iSamp],"l");

    if (iSamp==0) {
      hJetResPt[iSamp]->GetXaxis()->SetTitle(xlabel);
      hJetResPt[iSamp]->GetYaxis()->SetTitle(ylabel);
      hJetResPt[iSamp]->SetTitle("");
      //hJetResPt[iSamp]->GetYaxis()->SetRangeUser(-0.5,3);
      hJetResPt[iSamp]->Draw();

    }
    else {
      hJetResPt[iSamp]->Draw("same");
    }
  }
  //leg->Draw();
  c->SaveAs(TString(pname)+TString(".png"));

  hPt[0]->GetXaxis()->SetTitle("p_{T}");
  hPt[0]->SetTitle("");
  hPt[0]->SetLineColor(sampleColors[0]);
  hPt[0]->SetMarkerColor(sampleColors[0]);
  hPt[0]->Draw("hist");
  c->SaveAs("pt.png");

  sprintf(xlabel, "generator level tau Eta");
  sprintf(ylabel, "(reco P_{T} - gen P_{T})/gen P_{T}");
  sprintf(pname, "tauJetResEta");
  leg->Clear();

  for (UInt_t iSamp=0; iSamp<nSamples; iSamp++) {

    hJetResEta[iSamp]->SetLineColor(sampleColors[iSamp]);
    hJetResEta[iSamp]->SetMarkerColor(sampleColors[iSamp]);
    hJetResEta[iSamp]->SetMarkerSize(1);
    leg->AddEntry(hJetResEta[iSamp], sampleTitles[iSamp],"l");

    if (iSamp==0) {
      hJetResEta[iSamp]->GetXaxis()->SetTitle(xlabel);
      hJetResEta[iSamp]->GetYaxis()->SetTitle(ylabel);
      hJetResEta[iSamp]->SetTitle("");
      //hJetResEta[iSamp]->GetYaxis()->SetRangeUser(-0.3,0.3);
      hJetResEta[iSamp]->Draw();
    }

    else {
      hJetResEta[iSamp]->Draw("same");
    }
  }

  //leg->Draw();
  c->SaveAs(TString(pname)+TString(".png"));
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
