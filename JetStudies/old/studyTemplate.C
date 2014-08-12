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

void studyTemplate(const TString conf="sample.conf") {

  const Int_t nSamples=3;

  vector<TString> sampleNames;
  vector<TString> sampleTitles;
  vector<Int_t> sampleColors;

  confParse(conf, sampleNames, sampleTitles, sampleColors);

  TH1D *hSample[nSamples];

  char hname[100];

  for(UInt_t iSamp=0; iSamp<nSamples; iSamp++) {

    sprintf(hname, "hSample_%s", sampleTitles[iSamp].Data()); hSample[iSamp]= new TH1D(hname, "", 4, 0, 4);

  }

  UInt_t eventNum;
  UInt_t bTag1, bTag2;
  UInt_t tauDecayCat1, tauDecayCat2;
  LorentzVector *genB1=0, *genB2=0, *recoB1=0, *recoB2=0;
  LorentzVector *genTau1=0, *genTau2=0, *genDecayTau1=0, *genDecayTau2=0, *recoTau1=0, *recoTau2=0;

  TFile *infile;
  TTree *intree;

  for (UInt_t iSamp=0; iSamp<nSamples; iSamp++) { // sample loop

    TString infilename = sampleNames[iSamp];
    cout << "Processing  " << infilename << " ..." << endl;
    infile = new TFile(infilename); assert(infile);
    intree = (TTree*) infile->Get("Events"); assert(intree);
 
    intree->SetBranchAddress("eventNum",       &eventNum);
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

    for(UInt_t iEntry=0; iEntry<intree->GetEntries(); iEntry++) { // entry loop
      intree->GetEntry(iEntry);

      hSample[iSamp]->Fill(tauDecayCat1);

    } // end entry loop
    delete infile;
    infile=0, intree=0;

  } // end sample loop

  char pname[100];
  char xlabel[100];
  char ylabel[100];
  
  TCanvas *c = new TCanvas("c", "c", 800, 600);
  TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
  leg->SetShadowColor(0);
  leg->SetFillColor(0);

  for (UInt_t iSamp=0; iSamp<nSamples; iSamp++) {

    sprintf(xlabel, "Tau Decay Cat");
    sprintf(ylabel, "Events");
    sprintf(pname, "taudecay_%s", sampleTitles[iSamp].Data());
    hSample[iSamp]->SetLineColor(sampleColors[iSamp]);
    leg->AddEntry(hSample[iSamp], sampleTitles[iSamp],"l");

    if (iSamp==0) {
      hSample[iSamp]->GetXaxis()->SetTitle(xlabel);
      hSample[iSamp]->GetYaxis()->SetTitle(ylabel);
      hSample[iSamp]->Draw();
    }
    else {
      if (hSample[iSamp]->GetMaximum() > hSample[0]->GetMaximum()) hSample[0]->GetYaxis()->SetRangeUser(0,hSample[iSamp]->GetMaximum()*1.2);
      hSample[iSamp]->Draw("same");
    }
  }
  leg->Draw();
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
