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
#include <iomanip>
#include <fstream>
#include "Math/LorentzVector.h"

#include "MitStyleRemix.hh"
#include "bJetScaleCorr.h"
#endif

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;

void tauStudy(const TString inputfile = "/afs/cern.ch/work/k/klawhorn/HHToBBTT/Conf4_10GeV/HHToBBTT_10GeVJets_14TeV.root") {

  // set up input tree

  LorentzVector *genTau1=0, *genTau2=0, *genJetTau1=0, *genJetTau2=0, *jetTau1=0, *jetTau2=0;

  Float_t genTau1pt, genTau2pt, genJetTau1pt, genJetTau2pt, jetTau1pt, jetTau2pt, corrJetTau1pt, corrJetTau2pt;
  Float_t genTau1eta, genTau2eta, genJetTau1eta, genJetTau2eta, jetTau1eta, jetTau2eta;

  TFile *infile = new TFile(inputfile, "READ"); assert(infile);
  TTree *intree = (TTree*) infile->Get("Events"); assert(intree);

  intree->SetBranchAddress("genTau1",        &genTau1);
  intree->SetBranchAddress("genTau2",        &genTau2);
  intree->SetBranchAddress("genJetTau1",     &genJetTau1);
  intree->SetBranchAddress("genJetTau2",     &genJetTau2);
  intree->SetBranchAddress("jetTau1",        &jetTau1);
  intree->SetBranchAddress("jetTau2",        &jetTau2);

  const Float_t PT_MAX = 500;
  const Float_t PT_MIN = 0;
  const Int_t PT_BINS = 50;
  const Float_t ETA_MAX = 2.4;
  const Float_t ETA_MIN = -2.4;
  const Int_t ETA_BINS = 16;
  
  TProfile *tauTagEffWithPT = new TProfile("tauTagEffWithPT", "Tau Tagging Efficiency", PT_BINS, PT_MIN, PT_MAX);
  TProfile *tauTagEffWithEta = new TProfile("tauTagEffWithEta", "Tau Tagging Efficiency", ETA_BINS, ETA_MIN, ETA_MAX);

  //TProfile *ptScaleWithPT = new TProfile("ptScaleWithPT", "P_{T} Scale Factor between gen Tau and reco Jet", PT_BINS, PT_MIN, PT_MAX);
  //TProfile *ptJetScaleWithPT = new TProfile("ptJetScaleWithPT", "P_{T} Scale Factor between gen Jet and reco Jet", PT_BINS, PT_MIN, PT_MAX);

  //TProfile *ptScaleWithEta = new TProfile("ptScaleWithEta", "P_{T} Scale Factor between gen Tau and reco Jet", ETA_BINS, ETA_MIN, ETA_MAX);
  //TProfile *ptJetScaleWithEta = new TProfile("ptJetScaleWithEta", "P_{T} Scale Factor between gen Jet and reco Jet", ETA_BINS, ETA_MIN, ETA_MAX);

  //TProfile *ptResWithPT = new TProfile("ptResWithPT", "P_{T} Resolution between gen Tau and reco Jet", PT_BINS, PT_MIN, PT_MAX);
  TProfile *ptJetResWithPT = new TProfile("ptJetResWithPT", "P_{T} Resolution: (genJet-recoJet)/(genJet)", PT_BINS, PT_MIN, PT_MAX);
  TProfile *ptJetResCorrWithPT = new TProfile("ptJetResCorrWithPT", "P_{T} Resolution: (genJet-recoJet)/(genJet)", PT_BINS, PT_MIN, PT_MAX);

  //TProfile *ptResWithEta = new TProfile("ptResWithEta", "P_{T} Resolution between gen Tau and reco Jet", ETA_BINS, ETA_MIN, ETA_MAX);
  TProfile *ptJetResWithEta = new TProfile("ptJetResWithEta", "P_{T} Resolution: (genJet-recoJet)/(genJet)", ETA_BINS, ETA_MIN, ETA_MAX);
  TProfile *ptJetResCorrWithEta = new TProfile("ptJetResCorrWithEta", "P_{T} Resolution: (genJet-recoJet)/(genJet)", ETA_BINS, ETA_MIN, ETA_MAX);
  
  for (Int_t iEntry=0; iEntry<intree->GetEntries(); iEntry++) {
    intree->GetEntry(iEntry);

    genTau1pt = genTau1->Pt();       genTau1eta = genTau1->Eta();
    genTau2pt = genTau2->Pt();       genTau2eta = genTau2->Eta();
    genJetTau1pt = genJetTau1->Pt(); genJetTau1eta = genJetTau1->Eta();
    genJetTau2pt = genJetTau2->Pt(); genJetTau2eta = genJetTau2->Eta();
    jetTau1pt = jetTau1->Pt();       jetTau1eta = jetTau1->Eta();
    jetTau2pt = jetTau2->Pt();       jetTau2eta = jetTau2->Eta();

    // jet corrections
    //if ( (jetTau1pt != 999) ) cout << "eta " << jetTau1eta << " pt " << jetTau1pt << " scaleCorr " << getJetScaleFactor(jetTau1pt, jetTau1eta) << endl;
    corrJetTau1pt = jetTau1pt*getJetScaleFactor(jetTau1pt, jetTau1eta);
    corrJetTau2pt = jetTau2pt*getJetScaleFactor(jetTau2pt, jetTau2eta);

    if ( (genTau1pt != 999) && (jetTau1pt != 999) ) {
      tauTagEffWithPT->Fill(genTau1pt, 1);
      //ptScaleWithPT->Fill(genTau1pt, jetTau1pt/genTau1pt);
      //ptJetScaleWithPT->Fill(genJetTau1pt, jetTau1pt/genJetTau1pt);
      //ptResWithPT->Fill(genTau1pt, (genTau1pt-jetTau1pt)/genTau1pt);
      ptJetResWithPT->Fill(genJetTau1pt, (genJetTau1pt-jetTau1pt)/genJetTau1pt);
      ptJetResCorrWithPT->Fill(genJetTau1pt, (genJetTau1pt-corrJetTau1pt)/genJetTau1pt);

      tauTagEffWithEta->Fill(genTau1eta, 1);
      //ptScaleWithEta->Fill(genTau1eta, jetTau1pt/genTau1pt);
      //ptJetScaleWithEta->Fill(genJetTau1eta, jetTau1pt/genJetTau1pt);
      //ptResWithEta->Fill(genTau1eta, (genTau1pt-jetTau1pt)/genTau1pt);
      ptJetResWithEta->Fill(genJetTau1eta, (genJetTau1pt-jetTau1pt)/genJetTau1pt);
      ptJetResCorrWithEta->Fill(genJetTau1eta, (genJetTau1pt-corrJetTau1pt)/genJetTau1pt);

    }
    else if ( (genTau1pt != 999) && (jetTau1pt == 999) ) {
      tauTagEffWithPT->Fill(genTau1pt, 0);
      tauTagEffWithEta->Fill(genTau1eta, 0);
    }
    if ( (genTau2pt != 999) && (jetTau2pt != 999) ) {
      tauTagEffWithPT->Fill(genTau2pt, 1);
      //ptScaleWithPT->Fill(genTau2pt, jetTau2pt/genTau2pt);
      //ptJetScaleWithPT->Fill(genJetTau2pt, jetTau2pt/genJetTau2pt);
      //ptResWithPT->Fill(genTau2pt, (genTau2pt-jetTau2pt)/genTau2pt);
      ptJetResWithPT->Fill(genJetTau2pt, (genJetTau2pt-jetTau2pt)/genJetTau2pt);
      ptJetResCorrWithPT->Fill(genJetTau2pt, (genJetTau2pt-corrJetTau2pt)/genJetTau2pt);

      tauTagEffWithEta->Fill(genTau2eta, 1);
      //ptScaleWithEta->Fill(genTau2eta, jetTau2pt/genTau2pt);
      //ptJetScaleWithEta->Fill(genJetTau2eta, jetTau2pt/genJetTau2pt);
      //ptResWithEta->Fill(genTau2eta, (genTau2pt-jetTau2pt)/genTau2pt);
      ptJetResWithEta->Fill(genJetTau2eta, (genJetTau2pt-jetTau2pt)/genJetTau2pt);
      ptJetResCorrWithEta->Fill(genJetTau2eta, (genJetTau2pt-corrJetTau2pt)/genJetTau2pt);

    }
    else if ( (genTau2pt != 999) && (jetTau2pt == 999) ) {
      tauTagEffWithPT->Fill(genTau2pt, 0);
      tauTagEffWithEta->Fill(genTau2eta, 0);
    }

  }

  TCanvas *c1 = MakeCanvas("c1", "Tau Tagging Efficiency", 800, 600);
  tauTagEffWithPT->SetMarkerStyle(1);
  tauTagEffWithPT->GetXaxis()->SetTitle("P_{T} of Generator Tau");
  tauTagEffWithPT->GetYaxis()->SetTitle("Efficiency");
  tauTagEffWithPT->Draw();
  
  c1->SaveAs("tauTagEffWithPT.png");

  TCanvas *c2 = MakeCanvas("c2", "Tau Tagging Efficiency", 800, 600);
  tauTagEffWithEta->SetMarkerStyle(1);
  tauTagEffWithEta->GetXaxis()->SetTitle("Eta of Generator Tau");
  tauTagEffWithEta->GetYaxis()->SetTitle("Efficiency");
  tauTagEffWithEta->Draw();

  c2->SaveAs("tauTagEffWithEta.png");
  /*  
  TCanvas *c3 = MakeCanvas("c3", "Tau P_{T} Scale (reco Tau Jet P_{T}/gen Tau P_{T})", 800, 600);
  ptScaleWithPT->SetMarkerStyle(1);
  ptScaleWithPT->GetXaxis()->SetTitle("P_{T} of Generator Tau");
  ptScaleWithPT->GetYaxis()->SetTitle("P_{T} Scale");
  ptScaleWithPT->Draw();

  c3->SaveAs("ptScaleWithPT.png");

  TCanvas *c4 = MakeCanvas("c4", "Tau P_{T} Scale (reco Tau Jet P_{T}/gen Tau P_{T})", 800, 600);
  ptScaleWithEta->SetMarkerStyle(1);
  ptScaleWithEta->GetXaxis()->SetTitle("Eta of Generator Tau");
  ptScaleWithEta->GetYaxis()->SetTitle("P_{T} Scale");
  ptScaleWithEta->Draw();

  c4->SaveAs("ptScaleWithEta.png");

  TCanvas *c5 = MakeCanvas("c5", "Tau P_{T} Scale (reco Tau Jet P_{T}/gen Tau Jet P_{T})", 800, 600);
  ptJetScaleWithPT->SetMarkerStyle(1);
  ptJetScaleWithPT->GetXaxis()->SetTitle("P_{T} of Generator Tau Jet");
  ptJetScaleWithPT->GetYaxis()->SetTitle("P_{T} Scale");
  ptJetScaleWithPT->Draw();

  c5->SaveAs("ptJetScaleWithPT.png");

  TCanvas *c6 = MakeCanvas("c6", "Tau P_{T} Scale (reco Tau Jet P_{T}/gen Tau Jet P_{T})", 800, 600);
  ptJetScaleWithEta->SetMarkerStyle(1);
  ptJetScaleWithEta->GetXaxis()->SetTitle("Eta of Generator Tau Jet");
  ptJetScaleWithEta->GetYaxis()->SetTitle("P_{T} Scale");
  ptJetScaleWithEta->Draw();

  c6->SaveAs("ptJetScaleWithEta.png");

  TCanvas *c7 = MakeCanvas("c7", "Tau P_{T} Resolution ((reco Tau Jet P_{T} - gen Tau P_{T})/gen Tau P_{T})", 800, 600);
  ptResWithPT->SetMarkerStyle(1);
  ptResWithPT->GetXaxis()->SetTitle("P_{T} of Generator Tau");
  ptResWithPT->GetYaxis()->SetTitle("P_{T} Resolution");
  ptResWithPT->Draw();

  c7->SaveAs("ptResWithPT.png");

  TCanvas *c8 = MakeCanvas("c8", "Tau P_{T} Resolution", 800, 600);
  ptResWithEta->SetMarkerStyle(1);
  ptResWithEta->GetXaxis()->SetTitle("Eta of Generator Tau");
  ptResWithEta->GetYaxis()->SetTitle("P_{T} Resolution");
  ptResWithEta->Draw();
  
  c8->SaveAs("ptResWithEta.png");

  TCanvas *c9 = MakeCanvas("c9", "Tau P_{T} Resolution", 800, 600);
  ptJetResWithPT->SetMarkerStyle(1);
  ptJetResWithPT->GetXaxis()->SetTitle("P_{T} of Generator Tau Jet");
  ptJetResWithPT->GetYaxis()->SetTitle("P_{T} Resolution");
  ptJetResWithPT->Draw();
  ptJetResCorrWithPT->SetMarkerStyle(1);
  ptJetResCorrWithPT->SetLineColor(2);
  ptJetResCorrWithPT->Draw("same");

  c9->SaveAs("ptJetResWithPT.png");

  TCanvas *c10 = MakeCanvas("c10", "Tau P_{T} Resolution ((reco Tau Jet P_{T} - gen Tau Jet P_{T})/gen Tau Jet P_{T})", 800, 600);
  ptJetResWithEta->SetMarkerStyle(1);
  ptJetResWithEta->GetXaxis()->SetTitle("Eta of Generator Tau Jet");
  ptJetResWithEta->GetYaxis()->SetTitle("P_{T} Resolution");
  ptJetResWithEta->Draw();
  ptJetResCorrWithEta->SetMarkerStyle(1);
  ptJetResCorrWithEta->SetLineColor(2);
  ptJetResCorrWithEta->Draw("same");

  c10->SaveAs("ptJetResWithEta.png");
  */  
}
