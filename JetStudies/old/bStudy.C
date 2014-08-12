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

void bStudy(const TString conf="sample.conf") {

  // set up and parse input information
  const Int_t nSamples=2;

  vector<TString> sampleNames;
  vector<TString> sampleTitles;
  vector<Int_t> sampleColors;

  confParse(conf, sampleNames, sampleTitles, sampleColors);

  TProfile *hEff1Pt[nSamples], *hEff2Pt[nSamples];
  TProfile *hEff1Eta[nSamples], *hEff2Eta[nSamples];

  TH1D *hJetResUnCorPt[nSamples];//, *hJetResCorrPt[nSamples];
  TH1D *hJetResUnCorEta[nSamples];//, *hJetResCorrEta[nSamples];

  char hname[100];

  // define kinematic/plotting constants
  const Float_t PT_MAX  = 500;
  const Float_t PT_MIN  = 0;
  const Int_t   PT_BIN  = 10;
  const Float_t ETA_MAX = 4.0;
  const Float_t ETA_MIN = -4.0;
  const Int_t   ETA_BIN = 10;

  //Double_t jetCorr1Pt, jetCorr2Pt;

  // define plotting objects
  for(UInt_t iSamp=0; iSamp<nSamples; iSamp++) {

    // reco efficiency with pt and eta
    sprintf(hname, "hEff1Pt_%s", sampleTitles[iSamp].Data()); hEff1Pt[iSamp]= new TProfile(hname, hname, PT_BIN, PT_MIN, PT_MAX);
    sprintf(hname, "hEff2Pt_%s", sampleTitles[iSamp].Data()); hEff2Pt[iSamp]= new TProfile(hname, hname, PT_BIN, PT_MIN, PT_MAX);
    sprintf(hname, "hEff1Eta_%s", sampleTitles[iSamp].Data()); hEff1Eta[iSamp]= new TProfile(hname, hname, ETA_BIN, ETA_MIN, ETA_MAX);
    sprintf(hname, "hEff2Eta_%s", sampleTitles[iSamp].Data()); hEff2Eta[iSamp]= new TProfile(hname, hname, ETA_BIN, ETA_MIN, ETA_MAX);

    sprintf(hname, "hJetResUnCorPt_%s", sampleTitles[iSamp].Data()); hJetResUnCorPt[iSamp] = new TProfile(hname, hname, PT_BIN, PT_MIN, PT_MAX);
    //sprintf(hname, "hJetResCorrPt_%s", sampleTitles[iSamp].Data()); hJetResCorrPt[iSamp] = new TProfile(hname, hname, PT_BIN, PT_MIN, PT_MAX);
    sprintf(hname, "hJetResUnCorEta_%s", sampleTitles[iSamp].Data()); hJetResUnCorEta[iSamp] = new TProfile(hname, hname, ETA_BIN, ETA_MIN, ETA_MAX);
    //sprintf(hname, "hJetResCorrEta_%s", sampleTitles[iSamp].Data()); hJetResCorrEta[iSamp] = new TProfile(hname, hname, ETA_BIN, ETA_MIN, ETA_MAX);

  }

  // Define variables to read in files
  UInt_t eventNum;
  UInt_t bTag1, bTag2;
  UInt_t tauDecayCat1, tauDecayCat2;
  LorentzVector *genB1=0, *genB2=0, *recoB1=0, *recoB2=0;
  LorentzVector *genTau1=0, *genTau2=0, *genDecayTau1=0, *genDecayTau2=0, *recoTau1=0, *recoTau2=0;

  TFile *infile;
  TTree *intree;

  for (UInt_t iSamp=0; iSamp<nSamples; iSamp++) { // sample loop

    TString infilename = sampleNames[iSamp];
    cout << "Processing  " << infilename << " " << sampleTitles[iSamp] << " ..." << endl;
    infile = new TFile(infilename); assert(infile);
    intree = (TTree*) infile->Get("Events"); assert(intree);
 
    //intree->SetBranchAddress("eventNum",       &eventNum);
    intree->SetBranchAddress("bTag1",          &bTag1);
    intree->SetBranchAddress("bTag2",          &bTag2);
    intree->SetBranchAddress("sGenB1",          &genB1);
    intree->SetBranchAddress("sGenB2",          &genB2);
    intree->SetBranchAddress("sRecoB1",         &recoB1);
    intree->SetBranchAddress("sRecoB2",         &recoB2);
    /*intree->SetBranchAddress("tauCat1",        &tauDecayCat1);
    intree->SetBranchAddress("tauCat2",        &tauDecayCat2);
    intree->SetBranchAddress("genTau1",        &genTau1);
    intree->SetBranchAddress("genTau2",        &genTau2);
    intree->SetBranchAddress("genDecayTau1",   &genDecayTau1);
    intree->SetBranchAddress("genDecayTau2",   &genDecayTau2);
    intree->SetBranchAddress("recoTau1",       &recoTau1);
    intree->SetBranchAddress("recoTau2",       &recoTau2);*/

    for(UInt_t iEntry=0; iEntry<intree->GetEntries(); iEntry++) { // entry loop
      intree->GetEntry(iEntry);

      // jet resolution

      if (recoB1->Pt()!=999) {

	//jetCorr1Pt=recoB1->Pt()*getJetScaleFactor(recoB1->Pt(), recoB1->Eta());
	//cout << recoB1->Pt() << ", " << genB1->Pt() << "; ";
	//jetCorr1Pt=recoB1->Pt();

	hJetResUnCorPt[iSamp]->Fill(genB1->Pt(),(recoB1->Pt()-genB1->Pt())/genB1->Pt());
	hJetResUnCorEta[iSamp]->Fill(genB1->Eta(),(recoB1->Pt()-genB1->Pt())/genB1->Pt());
	//hJetResCorrPt[iSamp]->Fill(genB1->Pt(),(jetCorr1Pt-genB1->Pt())/genB1->Pt());
	//hJetResCorrEta[iSamp]->Fill(genB1->Eta(),(jetCorr1Pt-genB1->Pt())/genB1->Pt());

      }

      if (recoB2->Pt()!=999) {

	//jetCorr2Pt=recoB2->Pt()*getJetScaleFactor(recoB2->Pt(), recoB2->Eta());
	//jetCorr2Pt=recoB2->Pt();

	hJetResUnCorPt[iSamp]->Fill(genB2->Pt(),(recoB2->Pt()-genB2->Pt())/genB2->Pt());
	hJetResUnCorEta[iSamp]->Fill(genB2->Eta(),(recoB2->Pt()-genB2->Pt())/genB2->Pt());
	//hJetResCorrPt[iSamp]->Fill(genB2->Pt(),(jetCorr2Pt-genB2->Pt())/genB2->Pt());
	//hJetResCorrEta[iSamp]->Fill(genB2->Eta(),(jetCorr2Pt-genB2->Pt())/genB2->Pt());

      }

      // tagging efficiency

      if ( bTag1>0) {
	hEff1Pt[iSamp]->Fill(genB1->Pt(),1);
	hEff2Pt[iSamp]->Fill(genB1->Pt(),0);
	hEff1Eta[iSamp]->Fill(genB1->Eta(),1);
	hEff2Eta[iSamp]->Fill(genB1->Eta(),0);

      }
      else if (bTag1 == 2) {
	hEff1Pt[iSamp]->Fill(genB1->Pt(),0);
	hEff2Pt[iSamp]->Fill(genB1->Pt(),1);
	hEff1Eta[iSamp]->Fill(genB1->Eta(),0);
	hEff2Eta[iSamp]->Fill(genB1->Eta(),1);

      }
      else if (bTag1 == 3) {
	hEff1Pt[iSamp]->Fill(genB1->Pt(),1);
	hEff2Pt[iSamp]->Fill(genB1->Pt(),1);
	hEff1Eta[iSamp]->Fill(genB1->Eta(),1);
	hEff2Eta[iSamp]->Fill(genB1->Eta(),1);

      }
      else {
	//cout << recoB1->Eta() << endl;
	hEff1Pt[iSamp]->Fill(genB1->Pt(),0);
	hEff2Pt[iSamp]->Fill(genB1->Pt(),0);
	hEff1Eta[iSamp]->Fill(genB1->Eta(),0);
	hEff2Eta[iSamp]->Fill(genB1->Eta(),0);

      }

      if ( bTag2>0) {
	hEff1Pt[iSamp]->Fill(genB2->Pt(),1);
	hEff2Pt[iSamp]->Fill(genB2->Pt(),0);
	hEff1Eta[iSamp]->Fill(genB2->Eta(),1);
	hEff2Eta[iSamp]->Fill(genB2->Eta(),0);

      }
      else if (bTag2 == 2) {
	hEff1Pt[iSamp]->Fill(genB2->Pt(),0);
	hEff2Pt[iSamp]->Fill(genB2->Pt(),1);
	hEff1Eta[iSamp]->Fill(genB2->Eta(),0);
	hEff2Eta[iSamp]->Fill(genB2->Eta(),1);

      }
      else if (bTag2 == 3) {
	hEff1Pt[iSamp]->Fill(genB2->Pt(),1);
	hEff2Pt[iSamp]->Fill(genB2->Pt(),1);
	hEff1Eta[iSamp]->Fill(genB2->Eta(),1);
	hEff2Eta[iSamp]->Fill(genB2->Eta(),1);

      }
      else {
	hEff1Pt[iSamp]->Fill(genB2->Pt(),0);
	hEff2Pt[iSamp]->Fill(genB2->Pt(),0);
	hEff1Eta[iSamp]->Fill(genB2->Eta(),0);
	hEff2Eta[iSamp]->Fill(genB2->Eta(),0);

      }

    } // end entry loop
    delete infile;
    infile=0, intree=0;

  } // end sample loop
  
  char pname[100];
  char xlabel[100];
  char ylabel[100];
  
  TCanvas *c = MakeCanvas("c", "c", 800, 600);
  TLegend *leg = new TLegend(0.4, 0.2, 0.6, 0.4);
  leg->SetShadowColor(0);
  leg->SetFillColor(0);
  /*
  // b jet resolution as a function of pt (uncorrected)
  sprintf(xlabel, "generator level b P_{T}");
  sprintf(ylabel, "(reco P_{T} - gen P_{T})/gen P_{T}");
  sprintf(pname, "bJetResUnCorPt");
  
  for (UInt_t iSamp=0; iSamp<nSamples; iSamp++) {

    hJetResUnCorPt[iSamp]->SetLineColor(sampleColors[iSamp]);
    hJetResUnCorPt[iSamp]->SetMarkerColor(sampleColors[iSamp]);
    hJetResUnCorPt[iSamp]->SetMarkerSize(1);
    leg->AddEntry(hJetResUnCorPt[iSamp], sampleTitles[iSamp],"l");

    if (iSamp==0) {
      hJetResUnCorPt[iSamp]->GetXaxis()->SetTitle(xlabel);
      hJetResUnCorPt[iSamp]->GetYaxis()->SetTitle(ylabel);
      hJetResUnCorPt[iSamp]->SetTitle("");
      //hJetResUnCorPt[iSamp]->GetYaxis()->SetRangeUser(-0.5,3);
      hJetResUnCorPt[iSamp]->Draw();
    }

    else {
      hJetResUnCorPt[iSamp]->Draw("same");
    }
  }

  //leg->Draw();
  c->SaveAs(TString(pname)+TString(".png"));*/
  /*
  // b jet resolution as a function of pt (corrected)
  sprintf(xlabel, "generator level b P_{T}");
  sprintf(ylabel, "(reco P_{T} - gen P_{T})/gen P_{T}");
  sprintf(pname, "bJetResCorrPt");
  leg->Clear();

  for (UInt_t iSamp=0; iSamp<nSamples; iSamp++) {

    hJetResCorrPt[iSamp]->SetLineColor(sampleColors[iSamp]);
    hJetResCorrPt[iSamp]->SetMarkerColor(sampleColors[iSamp]);
    hJetResCorrPt[iSamp]->SetMarkerSize(1);
    leg->AddEntry(hJetResCorrPt[iSamp], sampleTitles[iSamp],"l");

    if (iSamp==0) {
      hJetResCorrPt[iSamp]->GetXaxis()->SetTitle(xlabel);
      hJetResCorrPt[iSamp]->GetYaxis()->SetTitle(ylabel);
      hJetResCorrPt[iSamp]->SetTitle("");
      //hJetResCorrPt[iSamp]->GetYaxis()->SetRangeUser(-0.5,3);
      hJetResCorrPt[iSamp]->Draw();
    }

    else {
      hJetResCorrPt[iSamp]->Draw("same");
    }
  }
  */
  //leg->Draw();
  /*  c->SaveAs(TString(pname)+TString(".png"));


  // b jet resolution as a function of eta (uncorrected)
  sprintf(xlabel, "generator level b Eta");
  sprintf(ylabel, "(reco P_{T} - gen P_{T})/gen P_{T}");
  sprintf(pname, "bJetResUnCorEta");
  leg->Clear();
  
  for (UInt_t iSamp=0; iSamp<nSamples; iSamp++) {

    hJetResUnCorEta[iSamp]->SetLineColor(sampleColors[iSamp]);
    hJetResUnCorEta[iSamp]->SetMarkerColor(sampleColors[iSamp]);
    hJetResUnCorEta[iSamp]->SetMarkerSize(1);
    leg->AddEntry(hJetResUnCorEta[iSamp], sampleTitles[iSamp],"l");

    if (iSamp==0) {
      hJetResUnCorEta[iSamp]->GetXaxis()->SetTitle(xlabel);
      hJetResUnCorEta[iSamp]->GetYaxis()->SetTitle(ylabel);
      hJetResUnCorEta[iSamp]->SetTitle("");
      //hJetResUnCorEta[iSamp]->GetYaxis()->SetRangeUser(-0.3,0.3);
      hJetResUnCorEta[iSamp]->Draw();
    }

    else {
      hJetResUnCorEta[iSamp]->Draw("same");
    }
  }

  //leg->Draw();
  c->SaveAs(TString(pname)+TString(".png"));
  */
  //leg->SetX1NDC(0.45); leg->SetX2NDC(0.65);
  //leg->SetY1NDC(0.2); leg->SetY2NDC(0.4);
  /*
  // b jet resolution as a function of eta (corrected)
  sprintf(xlabel, "generator level b Eta");
  sprintf(ylabel, "(reco P_{T} - gen P_{T})/gen P_{T}");
  sprintf(pname, "bJetResCorrEta");
  leg->Clear();

  for (UInt_t iSamp=0; iSamp<nSamples; iSamp++) {

    hJetResCorrEta[iSamp]->SetLineColor(sampleColors[iSamp]);
    hJetResCorrEta[iSamp]->SetMarkerColor(sampleColors[iSamp]);
    hJetResCorrEta[iSamp]->SetMarkerSize(1);
    leg->AddEntry(hJetResCorrEta[iSamp], sampleTitles[iSamp],"l");

    if (iSamp==0) {
      hJetResCorrEta[iSamp]->GetXaxis()->SetTitle(xlabel);
      hJetResCorrEta[iSamp]->GetYaxis()->SetTitle(ylabel);
      hJetResCorrEta[iSamp]->SetTitle("");
      //hJetResCorrEta[iSamp]->GetYaxis()->SetRangeUser(-0.3,0.3);
      hJetResCorrEta[iSamp]->Draw();
    }

    else {
      hJetResCorrEta[iSamp]->Draw("same");
    }
  }

  //leg->Draw();
  c->SaveAs(TString(pname)+TString(".png"));
  */
  // b jet efficiency as function of pt ( tagger1 )
  sprintf(xlabel, "gen b P_{T}");
  sprintf(ylabel, "b-tag efficiency");
  sprintf(pname, "bTag1effPt");
  leg->Clear();
 
  for (UInt_t iSamp=0; iSamp<nSamples; iSamp++) {

    hEff1Pt[iSamp]->SetLineColor(sampleColors[iSamp]);
    hEff1Pt[iSamp]->SetMarkerColor(sampleColors[iSamp]);
    hEff1Pt[iSamp]->SetMarkerSize(1);
    leg->AddEntry(hEff1Pt[iSamp], sampleTitles[iSamp],"l");

    // cout << hEff1Pt[iSamp]->FindLastBinAbove(0,1) << endl;

    if (iSamp==0) {
      hEff1Pt[iSamp]->GetXaxis()->SetTitle(xlabel);
      hEff1Pt[iSamp]->GetYaxis()->SetTitle(ylabel);
      hEff1Pt[iSamp]->SetTitle("");
      //hEff1Pt[iSamp]->GetYaxis()->SetRangeUser(0,1.0);
      hEff1Pt[iSamp]->Draw();

    }
    else {
    //if (hEff1Pt[iSamp]->GetMaximum() > hEff1Pt[0]->GetMaximum()) hEff1Pt[0]->GetYaxis()->SetRangeUser(0,hEff1Pt[iSamp]->GetMaximum()*1.2);
      hEff1Pt[iSamp]->Draw("same");
    }
  }
  leg->Draw();
  c->SaveAs(TString(pname)+TString(".png"));
  /*
  // b jet efficiency as function of pt (tagger2)
  sprintf(xlabel, "generator level b P_{T}");
  sprintf(ylabel, "tag algorithm 2 efficiency");
  sprintf(pname, "bTag2effPt");
  leg->Clear();

  for (UInt_t iSamp=0; iSamp<nSamples; iSamp++) {

    hEff2Pt[iSamp]->SetLineColor(sampleColors[iSamp]);
    hEff2Pt[iSamp]->SetMarkerColor(sampleColors[iSamp]);
    hEff2Pt[iSamp]->SetMarkerSize(1);
    leg->AddEntry(hEff2Pt[iSamp], sampleTitles[iSamp],"l");

    if (iSamp==0) {
      hEff2Pt[iSamp]->GetXaxis()->SetTitle(xlabel);
      hEff2Pt[iSamp]->GetYaxis()->SetTitle(ylabel);
      hEff2Pt[iSamp]->SetTitle("");
      //hEff2Pt[iSamp]->GetYaxis()->SetRangeUser(0.0,1.0);
      hEff2Pt[iSamp]->Draw();
    }
    else {
      //if (hEff2Pt[iSamp]->GetMaximum() > hEff2Pt[0]->GetMaximum()) hEff2Pt[0]->GetYaxis()->SetRangeUser(0,hEff2Pt[iSamp]->GetMaximum()*1.2);
      hEff2Pt[iSamp]->Draw("same");
    }
  }
  leg->Draw();
  c->SaveAs(TString(pname)+TString(".png"));
  */
  // b jet efficiency as function of eta (tagger1)
  sprintf(xlabel, "gen b #eta");
  sprintf(ylabel, "b-tag efficiency");
  sprintf(pname, "bTag1effEta");
  leg->Clear();

  leg->SetX1NDC(0.4); leg->SetX2NDC(0.6);
  leg->SetY1NDC(0.2); leg->SetY2NDC(0.4);

  for (UInt_t iSamp=0; iSamp<nSamples; iSamp++) {

    hEff1Eta[iSamp]->SetLineColor(sampleColors[iSamp]);
    hEff1Eta[iSamp]->SetMarkerColor(sampleColors[iSamp]);
    hEff1Eta[iSamp]->SetMarkerSize(1);
    leg->AddEntry(hEff1Eta[iSamp], sampleTitles[iSamp],"l");

    if (iSamp==0) {
      hEff1Eta[iSamp]->GetXaxis()->SetTitle(xlabel);
      hEff1Eta[iSamp]->GetYaxis()->SetTitle(ylabel);
      hEff1Eta[iSamp]->SetTitle("");
      //hEff1Eta[iSamp]->GetYaxis()->SetRangeUser(0.3,0.5);
      hEff1Eta[iSamp]->Draw();
    }
    else {
      //if (hEff1Eta[iSamp]->GetMaximum() > hEff1Eta[0]->GetMaximum()) hEff1Eta[0]->GetYaxis()->SetRangeUser(0,hEff1Eta[iSamp]->GetMaximum()*1.2);
      hEff1Eta[iSamp]->Draw("same");
    }
  }
  leg->Draw();
  c->SaveAs(TString(pname)+TString(".png"));
  /*
  // b jet efficiency as function of eta (tagger2)
  sprintf(xlabel, "generator level b eta");
  sprintf(ylabel, "tag algorithm 2 efficiency");
  sprintf(pname, "bTag2effEta");
  leg->Clear();

  for (UInt_t iSamp=0; iSamp<nSamples; iSamp++) {

    hEff2Eta[iSamp]->SetLineColor(sampleColors[iSamp]);
    hEff2Eta[iSamp]->SetMarkerColor(sampleColors[iSamp]);
    hEff2Eta[iSamp]->SetMarkerSize(1);
    leg->AddEntry(hEff2Eta[iSamp], sampleTitles[iSamp],"l");

    if (iSamp==0) {
      hEff2Eta[iSamp]->GetXaxis()->SetTitle(xlabel);
      hEff2Eta[iSamp]->GetYaxis()->SetTitle(ylabel);
      hEff2Eta[iSamp]->SetTitle("");
      //hEff2Eta[iSamp]->GetYaxis()->SetRangeUser(0.3,0.5);
      hEff2Eta[iSamp]->Draw();
    }
    else {
      //if (hEff2Eta[iSamp]->GetMaximum() > hEff2Eta[0]->GetMaximum()) hEff2Eta[0]->GetYaxis()->SetRangeUser(0,hEff2Eta[iSamp]->GetMaximum()*1.2);
      hEff2Eta[iSamp]->Draw("same");
    }
  }
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
