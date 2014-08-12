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

void delphesB(const TString conf="new.conf") {

  // set up and parse input information
  const Int_t nSamples=2;

  vector<TString> sampleNames;
  vector<TString> sampleTitles;
  vector<Int_t> sampleColors;

  confParse(conf, sampleNames, sampleTitles, sampleColors);

  TProfile *hEff1PtB[nSamples], *hEff2PtB[nSamples];
  TProfile *hEff1PtE[nSamples], *hEff2PtE[nSamples];

  char hname[100];

  // define kinematic/plotting constants
  const Float_t PT_MAX  = 400;
  const Float_t PT_MIN  = 0;
  const Int_t   PT_BIN  = 8;
  const Float_t ETA_MAX = 4.0;
  const Float_t ETA_TRANS = 1.2;

  // define plotting objects
  for(UInt_t iSamp=0; iSamp<nSamples; iSamp++) {

    // reco efficiency with pt and eta
    sprintf(hname, "hEff1PtB_%s", sampleTitles[iSamp].Data()); hEff1PtB[iSamp]= new TProfile(hname, hname, PT_BIN, PT_MIN, PT_MAX);
    sprintf(hname, "hEff2PtB_%s", sampleTitles[iSamp].Data()); hEff2PtB[iSamp]= new TProfile(hname, hname, PT_BIN, PT_MIN, PT_MAX);
    sprintf(hname, "hEff1PtE_%s", sampleTitles[iSamp].Data()); hEff1PtE[iSamp]= new TProfile(hname, hname, PT_BIN, PT_MIN, PT_MAX);
    sprintf(hname, "hEff2PtE_%s", sampleTitles[iSamp].Data()); hEff2PtE[iSamp]= new TProfile(hname, hname, PT_BIN, PT_MIN, PT_MAX);
  }

  // Define variables to read in files
  UInt_t bTag1, bTag2;
  LorentzVector *genB1=0, *genB2=0, *recoB1=0, *recoB2=0;

  TFile *infile;
  TTree *intree;

  for (UInt_t iSamp=0; iSamp<nSamples; iSamp++) { // sample loop

    TString infilename = sampleNames[iSamp];
    cout << "Processing  " << infilename << " " << sampleTitles[iSamp] << " ..." << endl;
    infile = new TFile(infilename); assert(infile);
    intree = (TTree*) infile->Get("Events"); assert(intree);
 
    intree->SetBranchAddress("bTag1",           &bTag1);
    intree->SetBranchAddress("bTag2",           &bTag2);
    intree->SetBranchAddress("sGenB1",          &genB1);
    intree->SetBranchAddress("sGenB2",          &genB2);
    intree->SetBranchAddress("sRecoB1",         &recoB1);
    intree->SetBranchAddress("sRecoB2",         &recoB2);

    for(UInt_t iEntry=0; iEntry<intree->GetEntries(); iEntry++) { // entry loop
      intree->GetEntry(iEntry);

      if (fabs(genB1->Eta())<ETA_TRANS) {
	if ( bTag1 == 1 ) {
	  hEff1PtB[iSamp]->Fill(genB1->Pt(),1);
	  hEff2PtB[iSamp]->Fill(genB1->Pt(),0);
	}
	else if (bTag1 == 2) {
	  hEff1PtB[iSamp]->Fill(genB1->Pt(),0);
	  hEff2PtB[iSamp]->Fill(genB1->Pt(),1);
	}
	else if (bTag1 == 3) {
	  hEff1PtB[iSamp]->Fill(genB1->Pt(),1);
	  hEff2PtB[iSamp]->Fill(genB1->Pt(),1);
	  
	}
	else {
	  hEff1PtB[iSamp]->Fill(genB1->Pt(),0);
	  hEff2PtB[iSamp]->Fill(genB1->Pt(),0);
	}
      }
      else if (fabs(genB1->Eta())<ETA_MAX) {
	if ( bTag1==1 ) {
	  hEff1PtE[iSamp]->Fill(genB1->Pt(),1);
	  hEff2PtE[iSamp]->Fill(genB1->Pt(),0);
	}
	else if (bTag1 == 2) {
	  hEff1PtE[iSamp]->Fill(genB1->Pt(),0);
	  hEff2PtE[iSamp]->Fill(genB1->Pt(),1);
	}
	else if (bTag1 == 3) {
	  hEff1PtE[iSamp]->Fill(genB1->Pt(),1);
	  hEff2PtE[iSamp]->Fill(genB1->Pt(),1);
	  
	}
	else {
	  hEff1PtE[iSamp]->Fill(genB1->Pt(),0);
	  hEff2PtE[iSamp]->Fill(genB1->Pt(),0);
	}
      }
      
      if (fabs(genB2->Eta())<ETA_TRANS) {
	if ( bTag2 == 1) {
	  hEff1PtB[iSamp]->Fill(genB2->Pt(),1);
	  hEff2PtB[iSamp]->Fill(genB2->Pt(),0);	  
	}
	else if (bTag2 == 2) {
	  hEff1PtB[iSamp]->Fill(genB2->Pt(),0);
	  hEff2PtB[iSamp]->Fill(genB2->Pt(),1);
	}
	else if (bTag2 == 3) {
	  hEff1PtB[iSamp]->Fill(genB2->Pt(),1);
	  hEff2PtB[iSamp]->Fill(genB2->Pt(),1);
	}
	else {
	  hEff1PtB[iSamp]->Fill(genB2->Pt(),0);
	  hEff2PtB[iSamp]->Fill(genB2->Pt(),0);
	}
      }
      else if (fabs(genB2->Eta())<ETA_MAX) {
	if ( bTag2 == 1) {
	  hEff1PtE[iSamp]->Fill(genB2->Pt(),1);
	  hEff2PtE[iSamp]->Fill(genB2->Pt(),0);	  
	}
	else if (bTag2 == 2) {
	  hEff1PtE[iSamp]->Fill(genB2->Pt(),0);
	  hEff2PtE[iSamp]->Fill(genB2->Pt(),1);
	}
	else if (bTag2 == 3) {
	  hEff1PtE[iSamp]->Fill(genB2->Pt(),1);
	  hEff2PtE[iSamp]->Fill(genB2->Pt(),1);
	}
	else {
	  hEff1PtE[iSamp]->Fill(genB2->Pt(),0);
	  hEff2PtE[iSamp]->Fill(genB2->Pt(),0);
	}
      }

    } // end entry loop
    
    cout << hEff1PtE[iSamp]->GetEntries() << endl;
    cout << hEff2PtE[iSamp]->GetEntries() << endl;
    cout << hEff1PtB[iSamp]->GetEntries() << endl;
    cout << hEff2PtB[iSamp]->GetEntries() << endl;    

    delete infile;
    infile=0, intree=0;

  } // end sample loop
  /*
  char pname[100];
  char xlabel[100];
  char ylabel[100];
  
  TCanvas *c = MakeCanvas("c", "c", 800, 600);
  TLegend *leg = new TLegend(0.2, 0.2, 0.35, 0.35);
  leg->SetShadowColor(0);
  leg->SetFillColor(0);

  // b jet efficiency as function of pt ( tagger1 )
  sprintf(xlabel, "gen b P_{T}");
  sprintf(ylabel, "b-tag 1 eff (barrel)");
  sprintf(pname, "bTag1pt_Bar");
  leg->Clear();

  TF1 *effxn = new TF1("effxn", "[0]*tanh([1]*x+[2])", 15, 500);

  effxn->SetParameter(0, 0.629858);
  effxn->SetParameter(1, 0.0166188);
  effxn->SetParameter(2, 0.300119);

  for (UInt_t iSamp=0; iSamp<nSamples; iSamp++) {

    hEff1PtB[iSamp]->SetLineColor(sampleColors[iSamp]);
    hEff1PtB[iSamp]->SetMarkerColor(sampleColors[iSamp]);
    hEff1PtB[iSamp]->SetMarkerSize(1);
    leg->AddEntry(hEff1PtB[iSamp], sampleTitles[iSamp],"l");

    if (iSamp==0) {
      hEff1PtB[iSamp]->GetXaxis()->SetTitle(xlabel);
      hEff1PtB[iSamp]->GetYaxis()->SetTitle(ylabel);
      hEff1PtB[iSamp]->SetTitle("");
      hEff1PtB[iSamp]->GetYaxis()->SetRangeUser(0.0, 1.0);
      hEff1PtB[iSamp]->Draw();
    }
    else {
      hEff1PtB[iSamp]->Draw("same");
    }
    effxn->Draw("same");
  }
  leg->Draw();
  c->SaveAs(TString(pname)+TString(".png"));

  // b jet efficiency as function of pt (tagger2)
  sprintf(xlabel, "generator level b P_{T}");
  sprintf(ylabel, "b-tag 2 eff (barrel)");
  sprintf(pname, "bTag2pt_Bar");
  leg->Clear();

  effxn->SetParameter(0, 0.629858);
  effxn->SetParameter(1, 0.0166188);
  effxn->SetParameter(2, 0.300119);  

  for (UInt_t iSamp=0; iSamp<nSamples; iSamp++) {

    hEff2PtB[iSamp]->SetLineColor(sampleColors[iSamp]);
    hEff2PtB[iSamp]->SetMarkerColor(sampleColors[iSamp]);
    hEff2PtB[iSamp]->SetMarkerSize(1);
    leg->AddEntry(hEff2PtB[iSamp], sampleTitles[iSamp],"l");

    if (iSamp==0) {
      hEff2PtB[iSamp]->GetXaxis()->SetTitle(xlabel);
      hEff2PtB[iSamp]->GetYaxis()->SetTitle(ylabel);
      hEff2PtB[iSamp]->SetTitle("");
      hEff2PtB[iSamp]->GetYaxis()->SetRangeUser(0.0, 1.0);
      hEff2PtB[iSamp]->Draw();
    }
    else {
      hEff2PtB[iSamp]->Draw("same");
    }
  }
  effxn->Draw("same");

  leg->Draw();
  c->SaveAs(TString(pname)+TString(".png"));

  // b jet efficiency as function of pt ( tagger1 )
  sprintf(xlabel, "gen b P_{T}");
  sprintf(ylabel, "b-tag 1 eff (endcap)");
  sprintf(pname, "bTag1pt_End");
  leg->Clear();

  effxn->SetParameter(0, 0.584522);
  effxn->SetParameter(1, 0.0144387);
  effxn->SetParameter(2, 0.397034);

  for (UInt_t iSamp=0; iSamp<nSamples; iSamp++) {

    hEff1PtE[iSamp]->SetLineColor(sampleColors[iSamp]);
    hEff1PtE[iSamp]->SetMarkerColor(sampleColors[iSamp]);
    hEff1PtE[iSamp]->SetMarkerSize(1);
    leg->AddEntry(hEff1PtE[iSamp], sampleTitles[iSamp],"l");

    if (iSamp==0) {
      hEff1PtE[iSamp]->GetXaxis()->SetTitle(xlabel);
      hEff1PtE[iSamp]->GetYaxis()->SetTitle(ylabel);
      hEff1PtE[iSamp]->SetTitle("");
      hEff1PtE[iSamp]->GetYaxis()->SetRangeUser(0.0, 1.0);
      hEff1PtE[iSamp]->Draw();
    }
    else {
      hEff1PtE[iSamp]->Draw("same");
    }
    effxn->Draw("same");
  }
  leg->Draw();
  c->SaveAs(TString(pname)+TString(".png"));

  // b jet efficiency as function of pt (tagger2)
  sprintf(xlabel, "generator level b P_{T}");
  sprintf(ylabel, "b-tag 2 eff (endcap)");
  sprintf(pname, "bTag2pt_End");
  leg->Clear();

  for (UInt_t iSamp=0; iSamp<nSamples; iSamp++) {

    hEff2PtE[iSamp]->SetLineColor(sampleColors[iSamp]);
    hEff2PtE[iSamp]->SetMarkerColor(sampleColors[iSamp]);
    hEff2PtE[iSamp]->SetMarkerSize(1);
    leg->AddEntry(hEff2PtE[iSamp], sampleTitles[iSamp],"l");

    if (iSamp==0) {
      hEff2PtE[iSamp]->GetXaxis()->SetTitle(xlabel);
      hEff2PtE[iSamp]->GetYaxis()->SetTitle(ylabel);
      hEff2PtE[iSamp]->SetTitle("");
      hEff2PtE[iSamp]->GetYaxis()->SetRangeUser(0.0, 1.0);
      hEff2PtE[iSamp]->Draw();
    }
    else {
      hEff2PtE[iSamp]->Draw("same");
    }
  }
  effxn->Draw("same");

  leg->Draw();
  c->SaveAs(TString(pname)+TString(".png"));
  */
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
