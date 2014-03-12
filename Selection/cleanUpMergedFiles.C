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
  Int_t nEvents;
  Float_t eventWeight;

  Float_t met, metPhi;

  Int_t eventType;
  UInt_t tauCat1=0, tauCat2=0;
  UInt_t bTag1=0, bTag2=0;

  LorentzVector *sRecoTau1=0, *sRecoTau2=0;
  LorentzVector *sGenJetTau1=0, *sGenJetTau2=0;
  LorentzVector *sGenTau1=0, *sGenTau2=0;

  LorentzVector *sRecoB1=0, *sRecoB2=0;
  LorentzVector *sGenJetB1=0, *sGenJetB2=0;
  LorentzVector *sGenB1=0, *sGenB2=0;

  LorentzVector *sRecoExtra=0;

  TFile* infile = new TFile(infilename); assert(infile);
  TTree* intree = (TTree*) infile->Get("Events"); assert(intree);

  inTree->SetBranchAddress("eventWeight",    &eventWeight);   // event weight from cross-section and Event->Weight
  inTree->SetBranchAddress("eventType",      &eventType);     // event type (0=signal, 1=tt, 2=zh, 3=other)
  inTree->SetBranchAddress("tauCat1",        &tauCat1);       // leading tau final state - jet, muon, electron
  inTree->SetBranchAddress("tauCat2",        &tauCat2);       // second tau final state - jet, muon, electron
  inTree->SetBranchAddress("bTag1",          &bTag1);         // leading b-jet tag from delphes
  inTree->SetBranchAddress("bTag2",          &bTag2);         // second b-jet tag from delphes
  inTree->SetBranchAddress("met",            &met,);          // missing transverse energy
  inTree->SetBranchAddress("metPhi",         &metPhi);        // missing transverse energy phi
  inTree->SetBranchAddress("mt2",            &mt2);           // "stransverse mass"
  inTree->SetBranchAddress("sGenTau1",       &sGenTau1);      // 4-vector for generator leading tau
  inTree->SetBranchAddress("sGenTau2",       &sGenTau2);      // 4-vector for generator second tau
  inTree->SetBranchAddress("sGenB1",         &sGenB1);        // 4-vector for generator leading b-jet
  inTree->SetBranchAddress("sGenB2",         &sGenB2);        // 4-vector for generator second b-jet
  inTree->SetBranchAddress("sGenJetTau1",    &sGenJetTau1);   // 4-vector for generator leading tau
  inTree->SetBranchAddress("sGenJetTau2",    &sGenJetTau2);   // 4-vector for generator second tau
  inTree->SetBranchAddress("sGenJetB1",      &sGenJetB1);     // 4-vector for generator leading b-jet
  inTree->SetBranchAddress("sGenJetB2",      &sGenJetB2);     // 4-vector for generator second b-jet
  inTree->SetBranchAddress("sRecoTau1",      &sRecoTau1);     // 4-vector for reconstructed leading tau
  inTree->SetBranchAddress("sRecoTau2",      &sRecoTau2);     // 4-vector for reconstructed second tau
  inTree->SetBranchAddress("sRecoB1",        &sRecoB1);       // 4-vector for reconstructed leading b-jet
  inTree->SetBranchAddress("sRecoB2",        &sRecoB2);       // 4-vector for reconstructed second b-jet
  inTree->SetBranchAddress("sRecoExtra",     &sRecoExtra);    // 4-vector for reconstructed extra jet

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
  outTree->Branch("eventType",      &eventType,      "eventType/i");    // event type (0=signal, 1=tt, 2=zh, 3=other)
  outTree->Branch("tauCat1",        &tauCat1,        "tauCat1/i");      // leading tau final state - jet, muon, electron
  outTree->Branch("tauCat2",        &tauCat2,        "tauCat2/i");      // second tau final state - jet, muon, electron
  outTree->Branch("bTag1",          &bTag1,          "bTag1/i");        // leading b-jet tag from delphes
  outTree->Branch("bTag2",          &bTag2,          "bTag2/i");        // second b-jet tag from delphes
  outTree->Branch("met",            &met,            "met/f");          // missing transverse energy
  outTree->Branch("metPhi",         &metPhi,         "metPhi/f");       // missing transverse energy phi
  outTree->Branch("mt2",            &mt2,            "mt2/D");          // "stransverse mass"
  outTree->Branch("sGenTau1",       "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sGenTau1);      // 4-vector for generator leading tau
  outTree->Branch("sGenTau2",       "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sGenTau2);      // 4-vector for generator second tau
  outTree->Branch("sGenB1",         "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sGenB1);        // 4-vector for generator leading b-jet
  outTree->Branch("sGenB2",         "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sGenB2);        // 4-vector for generator second b-jet
  outTree->Branch("sGenJetTau1",    "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sGenJetTau1);   // 4-vector for generator leading tau
  outTree->Branch("sGenJetTau2",    "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sGenJetTau2);   // 4-vector for generator second tau
  outTree->Branch("sGenJetB1",      "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sGenJetB1);     // 4-vector for generator leading b-jet
  outTree->Branch("sGenJetB2",      "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sGenJetB2);     // 4-vector for generator second b-jet
  outTree->Branch("sRecoTau1",      "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sRecoTau1);     // 4-vector for reconstructed leading tau
  outTree->Branch("sRecoTau2",      "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sRecoTau2);     // 4-vector for reconstructed second tau
  outTree->Branch("sRecoB1",        "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sRecoB1);       // 4-vector for reconstructed leading b-jet
  outTree->Branch("sRecoB2",        "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sRecoB2);       // 4-vector for reconstructed second b-jet
  outTree->Branch("sRecoExtra",     "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sRecoExtra);    // 4-vector for reconstructed extra jet

  for(UInt_t iEntry=0; iEntry<intree->GetEntries(); iEntry++) { // entry loop
    intree->GetEntry(iEntry);

    eventWeight=eventWeight/Float_t(totalEvents);

    outTree->Fill();

  }
  outFile->Write();
  outFile->Close();


}
