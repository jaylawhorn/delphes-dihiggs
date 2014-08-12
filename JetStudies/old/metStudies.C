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
#include <sstream>
#include "Math/LorentzVector.h"
#include "TLorentzVector.h"
#include "modules/Delphes.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "classes/DelphesClasses.h"
#include "../Utils/CPlot.hh"
#include "../Utils/MitStyleRemix.hh"

#endif

void metStudies() {

  TString inlist="hgg.txt";

  TChain chain("Delphes");

  ifstream ifs;
  ifs.open(inlist.Data()); assert(ifs.is_open());
  string line;
  TString fname;
  while(getline(ifs,line)) {
    stringstream ss(line);
    ss >> fname;
    chain.Add(fname);
  }
  ifs.close();

  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  TClonesArray *branchMET =treeReader->UseBranch("MissingET");
  TClonesArray *branchGenMET =treeReader->UseBranch("GenMissingET");
  TClonesArray *branchPuppiMET =treeReader->UseBranch("PuppiMissingET");
  TClonesArray *branchPileupMET =treeReader->UseBranch("PileUpJetIDMissingET");

  MissingET *met=0;

  //TH1D *hMetGen = new TH1D("hMetGen", "hMetGen", 30, 0, 300);
  //TH1D *hMet = new TH1D("hMet", "hMet", 30, 0, 300);
  //TH1D *hMetPuppi = new TH1D("hMetPuppi", "hMetPuppi", 30, 0, 300);
  //TH1D *hMetPU = new TH1D("hMetPU", "hMetPU", 30, 0, 300);

  TH1D *hMetXGen = new TH1D("hMetXGen", "hMetXGen", 60, -300, 300);
  TH1D *hMetX = new TH1D("hMetX", "hMetX", 60, -300, 300);
  TH1D *hMetXPuppi = new TH1D("hMetXPuppi", "hMetXPuppi", 60, -300, 300);
  TH1D *hMetXPU = new TH1D("hMetXPU", "hMetXPU", 60, -300, 300);

  TH1D *hMetYGen = new TH1D("hMetYGen", "hMetYGen", 60, -300, 300);
  TH1D *hMetY = new TH1D("hMetY", "hMetY", 60, -300, 300);
  TH1D *hMetYPuppi = new TH1D("hMetYPuppi", "hMetYPuppi", 60, -300, 300);
  TH1D *hMetYPU = new TH1D("hMetYPU", "hMetYPU", 60, -300, 300);

  /*TH1D *hRecoRes = new TH1D("hRecoRes", "hRecoRes", 28, -100, 600);
  TH1D *hPuppiRes = new TH1D("hPuppies", "hRecoRes", 28, -100, 600);
  TH1D *hPuRes = new TH1D("hPuRes", "hRecoRes", 28, -100, 600);

  TProfile *pReco = new TProfile("pReco", "pReco", 20, 0, 100);
  TProfile *pPuppi = new TProfile("pPuppi", "pPuppi", 20, 0, 100);
  TProfile *pPileup = new TProfile("pPileup", "pPileup", 20, 0, 100);*/

  Float_t genMetX=0, recoMetX=0, puppiMetX=0, puMetX=0;
  Float_t genMetY=0, recoMetY=0, puppiMetY=0, puMetY=0;
  
  for (Int_t iEntry=0; iEntry<numberOfEntries; iEntry++) { // entry loop
    treeReader->ReadEntry(iEntry);

    genMetX=0; recoMetX=0; puppiMetX=0; puMetX=0;
    genMetY=0; recoMetY=0; puppiMetY=0; puMetY=0;

    met = (MissingET*) branchGenMET->At(0); genMetX=met->MET*TMath::Cos(met->Phi); genMetY=met->MET*TMath::Sin(met->Phi); 
    met = (MissingET*) branchMET->At(0); recoMetX=met->MET*TMath::Cos(met->Phi); recoMetY=met->MET*TMath::Sin(met->Phi);
    met = (MissingET*) branchPuppiMET->At(0); puppiMetX=met->MET*TMath::Cos(met->Phi); puppiMetY=met->MET*TMath::Sin(met->Phi);
    met = (MissingET*) branchPileupMET->At(0); puMetX=met->MET*TMath::Cos(met->Phi); puMetY=met->MET*TMath::Sin(met->Phi);

    hMetXGen->Fill(genMetX);
    hMetX->Fill(recoMetX);
    hMetXPuppi->Fill(puppiMetX);
    hMetXPU->Fill(puMetX);

    hMetYGen->Fill(genMetY);
    hMetY->Fill(recoMetY);
    hMetYPuppi->Fill(puppiMetY);
    hMetYPU->Fill(puMetY);

    //hMetGen->Fill(genMet);
    //hMet->Fill(recoMet);
    //hMetPuppi->Fill(puppiMet);
    //hMetPU->Fill(puMet);

    /*   hRecoRes->Fill( (recoMet*TMath::Cos(-genMet)/genMet );
    hPuppiRes->Fill( (puppiMet-genMet)/genMet );
    hPuRes->Fill( (puMet-genMet)/genMet );

    pReco->Fill(genMet, recoMet/genMet);
    pPuppi->Fill(genMet, puppiMet/genMet);
    pPileup->Fill(genMet, puMet/genMet);*/

  }

  TCanvas *c1 = MakeCanvas("c1", "c1", 800, 800);

  TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
  leg->SetFillColor(0);
  
  hMetX->SetLineColor(kRed);
  hMetXPuppi->SetLineColor(kBlue);
  hMetXPU->SetLineColor(kGreen);
  hMetY->SetLineColor(kRed);
  hMetYPuppi->SetLineColor(kBlue);
  hMetYPU->SetLineColor(kGreen);
  leg->AddEntry(hMetXGen, "gen met", "l");
  leg->AddEntry(hMetX, "reco met", "l");
  leg->AddEntry(hMetXPuppi, "puppi met", "l");
  leg->AddEntry(hMetXPU, "pileup met", "l");

  hMetXGen->SetTitle("");
  hMetXGen->GetXaxis()->SetTitle("MET_{x}");
  hMetXGen->GetYaxis()->SetTitle("Events");
  hMetXGen->Draw();
  hMetX->Draw("same");
  hMetXPuppi->Draw("same");
  hMetXPU->Draw("same");
  leg->Draw("same");
  c1->SaveAs("met_x_raw.png");

  hMetYGen->SetTitle("");
  hMetYGen->GetXaxis()->SetTitle("MET_{y}");
  hMetYGen->GetYaxis()->SetTitle("Events");
  hMetYGen->Draw();
  hMetY->Draw("same");
  hMetYPuppi->Draw("same");
  hMetYPU->Draw("same");
  leg->Draw("same");
  c1->SaveAs("met_y_raw.png");

  /*
  leg->Clear();
  leg->AddEntry(hRecoRes, "reco met", "l");
  leg->AddEntry(hPuppiRes, "puppi met", "l");
  leg->AddEntry(hPuRes, "pileup met", "l");

  hRecoRes->SetLineColor(kRed);
  hPuppiRes->SetLineColor(kBlue);
  hPuRes->SetLineColor(kGreen);
  Float_t maxY=TMath::Max(TMath::Max(hRecoRes->GetMaximum(), hPuppiRes->GetMaximum()), hPuRes->GetMaximum());
  hRecoRes->SetTitle("");
  hRecoRes->GetXaxis()->SetTitle("(reco-gen)/gen");
  hRecoRes->GetYaxis()->SetTitle("Events");
  hRecoRes->GetYaxis()->SetRangeUser(0,1.1*maxY);
  hRecoRes->Draw();
  hPuppiRes->Draw("same");
  hPuRes->Draw("same");
  leg->Draw("same");
  c1->SaveAs("resolution.png");
    
  pReco->SetLineColor(kRed); pReco->SetMarkerColor(kRed);
  pPuppi->SetLineColor(kBlue); pPuppi->SetMarkerColor(kBlue);
  pPileup->SetLineColor(kGreen); pPileup->SetMarkerColor(kGreen);
  
  pReco->Draw();
  pPuppi->Draw("same");
  pPileup->Draw("same");
  leg->Draw("same");
  */
}
