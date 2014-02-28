//-------------------------------------------------------------------
// Clean up merged files from lxbtch
//
// execute with:
// root -l -q cleanUpMergedFiles(_infile_, _outfile_)
//
// Jay Lawhorn 11/4/13
//-------------------------------------------------------------------
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
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "Math/LorentzVector.h"

#endif

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;

void cleanUpMergedFiles(TString infilename="/afs/cern.ch/work/k/klawhorn/SnowmassSamples/PhaseII/Configuration4v2/Working/LL-4p-0-100-v1510_14TEV_temp.root",
			  TString outfilename="/afs/cern.ch/work/k/klawhorn/SnowmassSamples/PhaseII/Configuration4v2/Working/LL-4p-0-100-v1510_14TEV_clean.root") {

  // set up output variables and file
  UInt_t nEvents;
  UInt_t eventType;
  Float_t eventWeight;
  Float_t met, metPhi;
  Double_t mt2;
  UInt_t tauCat1=0, tauCat2=0;
  UInt_t bTag1=0, bTag2=0;
  LorentzVector *recoTau1=0, *recoTau2=0, *recoB1=0, *recoB2=0, *recoLeadJet=0, *recoExtraJet=0;
  LorentzVector *genTau1=0, *genTau2=0, *genDecayTau1=0, *genDecayTau2=0;
  //LorentzVector *genB1=0, *genB2=0, *recoB1=0, *recoB2=0, *boostJet=0, *genBoostJet=0;

  TFile* infile = new TFile(infilename); assert(infile);
  TTree* intree = (TTree*) infile->Get("Events"); assert(intree);

  intree->SetBranchAddress("eventWeight",    &eventWeight);
  intree->SetBranchAddress("eventType",      &eventType);
  intree->SetBranchAddress("tauCat1",        &tauCat1);
  intree->SetBranchAddress("tauCat2",        &tauCat2);
  intree->SetBranchAddress("bTag1",          &bTag1);
  intree->SetBranchAddress("bTag2",          &bTag2);
  intree->SetBranchAddress("met",            &met);
  intree->SetBranchAddress("metPhi",         &metPhi);
  intree->SetBranchAddress("mt2",            &mt2);
  intree->SetBranchAddress("recoTau1",       &recoTau1);     // 4-vector for reconstructed leading tau
  intree->SetBranchAddress("recoTau2",       &recoTau2);     // 4-vector for reconstructed second tau
  intree->SetBranchAddress("genTau1",        &genTau1);    
  intree->SetBranchAddress("genTau2",        &genTau2);    
  intree->SetBranchAddress("genDecayTau1",   &genDecayTau1);
  intree->SetBranchAddress("genDecayTau2",   &genDecayTau2);
  intree->SetBranchAddress("recoB1",         &recoB1);       // 4-vector for reconstructed leading b-jet
  intree->SetBranchAddress("recoB2",         &recoB2);       // 4-vector for reconstructed second b-jet
  intree->SetBranchAddress("recoLeadJet",    &recoLeadJet);  // 4-vector for reconstructed leading jet
  intree->SetBranchAddress("recoExtraJet",   &recoExtraJet); // 4-vector for reconstructed extra jet

  TTree* infotree = (TTree*) infile->Get("Info"); assert(infotree);
  infotree->SetBranchAddress("nEvents",      &nEvents);

  Int_t totalEvents=0;

  for (UInt_t iEntry=0; iEntry<infotree->GetEntries(); iEntry++) {
    infotree->GetEntry(iEntry);
    totalEvents+=nEvents;
  }

  TFile *outFile = new TFile(outfilename, "RECREATE");

  // tree to hold information about selected events
  TTree *outTree = new TTree("Events", "Events");
  outTree->Branch("eventWeight",    &eventWeight,    "eventWeight/f");  // event weight from cross-section and Event->Weight
  outTree->Branch("eventType",      &eventType,      "eventType/i");    // event type
  outTree->Branch("tauCat1",        &tauCat1,        "tauCat1/i");      // leading tau final state - jet, muon, electron
  outTree->Branch("tauCat2",        &tauCat2,        "tauCat2/i");      // second tau final state - jet, muon, electron
  outTree->Branch("bTag1",          &bTag1,          "bTag1/i");        // leading b-jet tag from delphes
  outTree->Branch("bTag2",          &bTag2,          "bTag2/i");        // second b-jet tag from delphes
  outTree->Branch("met",            &met,            "met/f");          // missing transverse energy
  outTree->Branch("metPhi",         &metPhi,         "metPhi/f");       // missing transverse energy phi
  outTree->Branch("mt2",            &mt2,            "mt2/D");          // "stransverse mass" variable
  outTree->Branch("recoTau1",       "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &recoTau1);     // 4-vector for reconstructed leading tau
  outTree->Branch("recoTau2",       "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &recoTau2);     // 4-vector for reconstructed second tau
  outTree->Branch("genTau1",        "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &genTau1);
  outTree->Branch("genTau2",        "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &genTau2);
  outTree->Branch("genDecayTau1",   "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &genDecayTau1);
  outTree->Branch("genDecayTau2",   "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &genDecayTau2);
  outTree->Branch("recoB1",         "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &recoB1);       // 4-vector for reconstructed leading b-jet
  outTree->Branch("recoB2",         "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &recoB2);       // 4-vector for reconstructed second b-jet
  outTree->Branch("recoLeadJet",    "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &recoLeadJet);  // 4-vector for reconstructed leading jet
  outTree->Branch("recoExtraJet",   "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &recoExtraJet); // 4-vector for reconstructed extra jet

  for(UInt_t iEntry=0; iEntry<intree->GetEntries(); iEntry++) { // entry loop
    intree->GetEntry(iEntry);

    eventWeight=eventWeight/Float_t(totalEvents);

    outTree->Fill();

  }
  outFile->Write();
  outFile->Close();


}
