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

Float_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 );
void confParse(const TString conf, vector<TString> &sampleNames, vector<TString> &sampleTitles, vector<Int_t> &sampleColors);

void genBstudy(const TString conf="gen.conf") {

  const Int_t nSamples=1;

  vector<TString> sampleNames;
  vector<TString> sampleTitles;
  vector<Int_t> sampleColors;

  confParse(conf, sampleNames, sampleTitles, sampleColors);

  TH1D *hBdR[nSamples];

  char hname[100];

  for(UInt_t iSamp=0; iSamp<nSamples; iSamp++) {

    sprintf(hname, "hBdR_%s", sampleTitles[iSamp].Data()); hBdR[iSamp]= new TH1D(hname, "", 32, 0, 8.0); hBdR[iSamp]->Sumw2();

  }

  Float_t eventWeight;
  LorentzVector *genB1=0, *genB2=0;
  LorentzVector *genTau1=0, *genTau2=0;

  TFile *infile;
  TTree *intree;

  for (UInt_t iSamp=0; iSamp<nSamples; iSamp++) { // sample loop

    TString infilename = sampleNames[iSamp];
    cout << "Processing  " << infilename << " ..." << endl;
    infile = new TFile(infilename); assert(infile);
    intree = (TTree*) infile->Get("Events"); assert(intree);
 
    intree->SetBranchAddress("eventWeight",    &eventWeight);
    intree->SetBranchAddress("genB1",          &genB1);
    intree->SetBranchAddress("genB2",          &genB2);
    intree->SetBranchAddress("genTau1",        &genTau1);
    intree->SetBranchAddress("genTau2",        &genTau2);

    for(UInt_t iEntry=0; iEntry<intree->GetEntries(); iEntry++) { // entry loop
      intree->GetEntry(iEntry);

      hBdR[iSamp]->Fill(deltaR(genB1->Eta(), genB2->Eta(), genB1->Phi(), genB2->Phi()), eventWeight);

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

    sprintf(xlabel, "#Deltar (bb)");
    sprintf(ylabel, "Events");
    sprintf(pname, "bDr_%s", sampleTitles[iSamp].Data());
    hBdR[iSamp]->SetLineColor(sampleColors[iSamp]);
    leg->AddEntry(hBdR[iSamp], sampleTitles[iSamp],"l");

    if (iSamp==0) {
      hBdR[iSamp]->SetLineColor(kBlack);
      hBdR[iSamp]->GetXaxis()->SetTitle(xlabel);
      hBdR[iSamp]->GetYaxis()->SetTitle(ylabel);
      Float_t scale= hBdR[iSamp]->Integral();
      hBdR[iSamp]->Scale(1/scale);
      hBdR[iSamp]->Draw("hist");
    }
    else {
      if (hBdR[iSamp]->GetMaximum() > hBdR[0]->GetMaximum()) hBdR[0]->GetYaxis()->SetRangeUser(0,hBdR[iSamp]->GetMaximum()*1.2);
      hBdR[iSamp]->Draw("same");
    }
  }
  //leg->Draw();

  c->SaveAs("bDr.png");
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
    cout << fname << title << color << endl;
    sampleNames.push_back(fname);
    sampleTitles.push_back(title);
    sampleColors.push_back(color);

  }
  ifs.close();

}

Float_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 ) {

  const Float_t pi = 3.14159265358979;

  Float_t etaDiff = (eta1-eta2);
  Float_t phiDiff = fabs(phi1-phi2);
  while ( phiDiff>pi ) phiDiff = fabs(phiDiff-2.0*pi);

  Float_t deltaRSquared = etaDiff*etaDiff + phiDiff*phiDiff;

  return TMath::Sqrt(deltaRSquared);

}
