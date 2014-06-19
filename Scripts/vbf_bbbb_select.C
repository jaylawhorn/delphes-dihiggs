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
#include "modules/Delphes.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "classes/DelphesClasses.h"
#include "mt2.hh"
#endif

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;

Float_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 );

Int_t puJetID( Float_t eta, Float_t meanSqDeltaR, Float_t beta );

void vbf_bbbb_select(const TString confname="test_vbf.txt",
	       const Float_t xsec=1.0,
	       const Float_t totalEvents=100,
	       const TString outputfile="test.root") {

  // declare constants
  const Int_t B_ID_CODE = 5;
  const Int_t H_ID_CODE = 25;

  const Float_t MAX_MATCH_DIST = 0.5;

 // read input input file
  TChain chain("Delphes");
  TString inputfile;
  
  ifstream ifs;
  ifs.open(confname);
  assert(ifs.is_open());
  string line;
  while (getline(ifs, line)) {
    stringstream ss(line);;
    ss >> inputfile;
    chain.Add(inputfile);    
  }
  ifs.close();
  
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchGenJet = treeReader->UseBranch("GenJet");
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");

  //set up loop variables
  GenParticle *genParticle=0;
  Jet *genJet=0;
  Jet *jet=0;

  // set up storage variables
  Jet         *jetB1=0, *jetB2=0, *jetB3=0, *jetB4=0;
  Jet         *jet1=0, *jet2=0;
  GenParticle *genB1=0, *genB2=0, *genB3=0, *genB4=0;
  GenParticle *genH1=0, *genH2=0;
  Jet *genJetB1=0, *genJetB2=0, *genJetB3=0, *genJetB4=0;
  Jet *genJetVBF1=0, *genJetVBF2=0;

  Int_t iB1=-1,       iB2=-1,       iB3=-1,       iB4=-1;
  Int_t iH1=-1,       iH2=-1;
  Int_t iJet1=-1,     iJet2=-1;
  Int_t iGenB1=-1,    iGenB2=-1,    iGenB3=-1,    iGenB4=-1;
  Int_t iGenJetB1=-1, iGenJetB2=-1, iGenJetB3=-1, iGenJetB4=-1;
  Int_t iGenJetVBF1=-1, iGenJetVBF2=-1;
  Int_t iHmatch1=-1,   iHmatch2=-1,  iHmatch3=-1,  iHmatch4=-1;

  LorentzVector *sRecoB1=0, *sRecoB2=0, *sRecoB3=0, *sRecoB4=0;
  LorentzVector *sGenJetB1=0, *sGenJetB2=0, *sGenJetB3=0, *sGenJetB4=0;
  LorentzVector *sGenJetVBF1=0, *sGenJetVBF2=0;
  LorentzVector *sGenB1=0, *sGenB2=0, *sGenB3=0, *sGenB4=0;
  LorentzVector *sGenH1=0, *sGenH2=0;
  LorentzVector *sRecoJet1=0, *sRecoJet2=0;

  TFile *outFile = new TFile(outputfile, "RECREATE");

  // tree to hold information about selected events
  TTree *outTree = new TTree("Events", "Events");

  outTree->Branch("iHmatch1",       &iHmatch1,                                                     "iHmatch1/i");   // which Higgs does b-jet 1 come from
  outTree->Branch("iHmatch2",       &iHmatch2,                                                     "iHmatch2/i");   // which Higgs does b-jet 2 come from
  outTree->Branch("iHmatch3",       &iHmatch3,                                                     "iHmatch3/i");   // which Higgs does b-jet 3 come from
  outTree->Branch("iHmatch4",       &iHmatch4,                                                     "iHmatch4/i");   // which Higgs does b-jet 4 come from

  outTree->Branch("sGenB1",         "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sGenB1);        // 4-vector for generator leading b-jet
  outTree->Branch("sGenB2",         "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sGenB2);        // 4-vector for generator b-jet
  outTree->Branch("sGenB3",         "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sGenB3);        // 4-vector for generator b-jet
  outTree->Branch("sGenB4",         "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sGenB4);        // 4-vector for generator b-jet

  outTree->Branch("sGenH1",         "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sGenH1);        // 4-vector for generator higgs
  outTree->Branch("sGenH2",         "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sGenH2);        // 4-vector for generator higgs

  outTree->Branch("sGenJetB1",      "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sGenJetB1);     // 4-vector for generator leading b-jet
  outTree->Branch("sGenJetB2",      "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sGenJetB2);     // 4-vector for generator b-jet
  outTree->Branch("sGenJetB3",      "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sGenJetB3);     // 4-vector for generator b-jet
  outTree->Branch("sGenJetB4",      "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sGenJetB4);     // 4-vector for generator b-jet

  outTree->Branch("sGenJetVBF1",    "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sGenJetVBF1);   // 4-vector for generator leading b-jet
  outTree->Branch("sGenJetVBF2",    "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sGenJetVBF2);   // 4-vector for generator b-jet

  outTree->Branch("sRecoB1",        "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sRecoB1);       // 4-vector for reconstructed leading b-jet
  outTree->Branch("sRecoB2",        "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sRecoB2);       // 4-vector for reconstructed b-jet
  outTree->Branch("sRecoB3",        "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sRecoB3);       // 4-vector for reconstructed b-jet
  outTree->Branch("sRecoB4",        "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sRecoB4);       // 4-vector for reconstructed b-jet

  outTree->Branch("sRecoJet1",      "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sRecoJet1);     // 4-vector for reconstructed leading forward-jet
  outTree->Branch("sRecoJet2",      "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sRecoJet2);     // 4-vector for reconstructed second forward-jet

  // define placeholder vector for things that don't exist
  LorentzVector nothing(-999,-999,0,-999);

  Int_t iMatched=0;
  Int_t iNot=0;
  Int_t iTwo=0;

  Int_t nM=0;

  for (Int_t iEntry=0; iEntry<numberOfEntries; iEntry++) { // entry loop
    treeReader->ReadEntry(iEntry);

    // ********************
    // RESET
    // ********************
    
    iB1=-1; iB2=-1; iB3=-1; iB4=-1;
    iGenB1=-1; iGenB2=-1; iGenB3=-1; iGenB4=-1;
    iGenJetB1=-1; iGenJetB2=-1; iGenJetB3=-1; iGenJetB4=-1;
    iGenJetVBF1=-1; iGenJetVBF2=-1;
    iH1=-1; iH2=-1;
    iJet1=-1; iJet2=-1;
    iHmatch1=-1; iHmatch2=-1; iHmatch3=-1; iHmatch4=-1;

    jet1=0; jet2=0;
    jetB1=0; jetB2=0; jetB3=0; jetB4=0;
    genB1=0; genB2=0; genB3=0; genB4=0;
    genJetB1=0; genJetB2=0; genJetB3=0; genJetB4=0;
    genJetVBF1=0; genJetVBF2=0;
    sGenB1=0; sGenB2=0; sGenB3=0; sGenB4=0;
    sGenJetB1=0; sGenJetB2=0; sGenJetB3=0; sGenJetB4=0;
    sGenJetVBF1=0; sGenJetVBF2=0;
    sGenH1=0; sGenH2=0;

    nM=0;

    for (Int_t iJet=0; iJet<branchJet->GetEntries(); iJet++) { // reconstructed jet loop
      jet = (Jet*) branchJet->At(iJet);

      if (fabs(jet->Eta)>4.0) continue;
      if (jet->PT<30) continue;
      if (puJetID(jet->Eta, jet->MeanSqDeltaR, jet->BetaStar)==1) continue;
      if (jet->BTag==0) continue;

      if ((jetB1)&&(deltaR(jet->Eta, jetB1->Eta, jet->Phi, jetB1->Phi) < MAX_MATCH_DIST)) continue;
      if ((jetB2)&&(deltaR(jet->Eta, jetB2->Eta, jet->Phi, jetB2->Phi) < MAX_MATCH_DIST)) continue;
      if ((jetB3)&&(deltaR(jet->Eta, jetB3->Eta, jet->Phi, jetB3->Phi) < MAX_MATCH_DIST)) continue;
      if ((jetB4)&&(deltaR(jet->Eta, jetB4->Eta, jet->Phi, jetB4->Phi) < MAX_MATCH_DIST)) continue;

      if (iB1==-1) {
	iB1=iJet; 
	jetB1 = (Jet*) branchJet->At(iB1); 
      }
      else if (jet->PT > jetB1->PT) {
	iB4=iB3;
	jetB4 = (Jet*) branchJet->At(iB4);
	iB3=iB2;
	jetB3 = (Jet*) branchJet->At(iB3);
	iB2=iB1; 
	jetB2 = (Jet*) branchJet->At(iB2); 
	iB1=iJet;
	jetB1 = (Jet*) branchJet->At(iB1);
      }
      else if (iB2==-1) { 
	iB2=iJet; 
	jetB2 = (Jet*) branchJet->At(iB2); 
      }
      else if (jet->PT > jetB2->PT) { 
	iB4=iB3;
        jetB4 = (Jet*) branchJet->At(iB4);
        iB3=iB2;
        jetB3 = (Jet*) branchJet->At(iB3);
	iB2=iJet; 
	jetB2 = (Jet*) branchJet->At(iB2); 
      }
      else if (iB3==-1) {
	iB3=iJet;
	jetB3 = (Jet*) branchJet->At(iB3);
      }
      else if (jet->PT > jetB3->PT) {
	iB4=iB3;
        jetB4 = (Jet*) branchJet->At(iB4);
        iB3=iB2;
        jetB3 = (Jet*) branchJet->At(iB3);
      }
      else if (iB4==-1) {
	iB4=iJet;
	jetB4 = (Jet*) branchJet->At(iB4);
      }
      else if (jet->PT > jetB4->PT) {
	iB4=iJet;
	jetB4 = (Jet*) branchJet->At(iB4);
      }

    } // end reco jet loop

    // get VBF jets
    for (Int_t iJet=0; iJet<branchJet->GetEntries(); iJet++) { // reconstructed jet loop
      jet = (Jet*) branchJet->At(iJet);

      if (fabs(jet->Eta)>4.7) continue;
      if (jet->PT<30) continue;

      if (puJetID(jet->Eta, jet->MeanSqDeltaR, jet->BetaStar)==1) continue;

      if ((jetB1)&&(deltaR(jet->Eta, jetB1->Eta, jet->Phi, jetB1->Phi) < MAX_MATCH_DIST)) continue;
      if ((jetB2)&&(deltaR(jet->Eta, jetB2->Eta, jet->Phi, jetB2->Phi) < MAX_MATCH_DIST)) continue;
      if ((jetB3)&&(deltaR(jet->Eta, jetB3->Eta, jet->Phi, jetB3->Phi) < MAX_MATCH_DIST)) continue;
      if ((jetB4)&&(deltaR(jet->Eta, jetB4->Eta, jet->Phi, jetB4->Phi) < MAX_MATCH_DIST)) continue;

      if (iJet1==-1) {
        iJet1=iJet;
        jet1 = (Jet*) branchJet->At(iJet1);
      }
      else if (jet->PT > jet1->PT) {
        iJet2=iJet1;
        jet2 = (Jet*) branchJet->At(iJet2);
        iJet1=iJet;
        jet1 = (Jet*) branchJet->At(iJet1);
      }
      else if (iJet2==-1) {
        iJet2=iJet;
        jet2 = (Jet*) branchJet->At(iJet2);
      }
      else if (jet->PT > jet2->PT) {
        iJet2=iJet;
        jet2 = (Jet*) branchJet->At(iJet2);
      }
    }

    if ( (!jetB1) || (!jetB2) || (!jetB3) || (!jetB4) || (!jet1) || (!jet2) ) continue;

    /*    cout << "stored jets" << endl;
    if (jetB1) cout << "1 " <<  jetB1->PT << " " << jetB1->Eta << endl;
    if (jetB2) cout << "2 " << jetB2->PT << " " << jetB2->Eta << endl;
    if (jetB3) cout << "3 " << jetB3->PT << " " << jetB3->Eta << endl;
    if (jetB4) cout << "4 " << jetB4->PT << " " << jetB4->Eta << endl;
    if (jet1) cout << "V1 " << jet1->PT << " " << jet1->Eta << endl;
    if (jet2) cout << "V2 " << jet2->PT << " " << jet2->Eta << endl;
    cout << endl;*/

    // fill 4-vector for leading b-jet
    LorentzVector vRecoB1(0,0,0,0);
    if (jetB1) {
      vRecoB1.SetPt(jetB1->PT);
      vRecoB1.SetEta(jetB1->Eta);
      vRecoB1.SetPhi(jetB1->Phi);
      vRecoB1.SetM(jetB1->Mass);
      sRecoB1 = &vRecoB1;
    }
    else sRecoB1 = &nothing;

    // fill 4-vector for b-jet
    LorentzVector vRecoB2(0,0,0,0);
    if (jetB2) {
      vRecoB2.SetPt(jetB2->PT);
      vRecoB2.SetEta(jetB2->Eta);
      vRecoB2.SetPhi(jetB2->Phi);
      vRecoB2.SetM(jetB2->Mass);
      sRecoB2 = &vRecoB2;
    }
    else sRecoB2 = &nothing;

    // fill 4-vector for b-jet
    LorentzVector vRecoB3(0,0,0,0);
    if (jetB3) {
      vRecoB3.SetPt(jetB3->PT);
      vRecoB3.SetEta(jetB3->Eta);
      vRecoB3.SetPhi(jetB3->Phi);
      vRecoB3.SetM(jetB3->Mass);
      sRecoB3 = &vRecoB3;
    }
    else sRecoB3 = &nothing;

    // fill 4-vector for b-jet
    LorentzVector vRecoB4(0,0,0,0);
    if (jetB4) {
      vRecoB4.SetPt(jetB4->PT);
      vRecoB4.SetEta(jetB4->Eta);
      vRecoB4.SetPhi(jetB4->Phi);
      vRecoB4.SetM(jetB4->Mass);
      sRecoB4 = &vRecoB4;
    }
    else sRecoB4 = &nothing;

    // ********************
    // GEN PARTICLES
    // ********************
    
    //cout << "B-QUARKS: " << endl;
    for (Int_t iParticle=0; iParticle<branchParticle->GetEntries(); iParticle++) { // generator particle loop
      genParticle = (GenParticle*) branchParticle->At(iParticle);

      if ( fabs(genParticle->PID) != B_ID_CODE ) continue;
      if (genParticle->Status != 3) continue;
      //cout << genParticle->PT << " " << genParticle->Eta << " " << genParticle->Phi << " " << genParticle->Status << endl;

      if ( deltaR(genParticle->Eta, vRecoB1.Eta(), genParticle->Phi, vRecoB1.Phi()) < MAX_MATCH_DIST ) {
	iGenB1 = iParticle;
	genB1 = (GenParticle*) branchParticle->At(iGenB1);
      }
      else if ( deltaR(genParticle->Eta, vRecoB2.Eta(), genParticle->Phi, vRecoB2.Phi()) < MAX_MATCH_DIST ) {
	iGenB2 = iParticle;
	genB2 = (GenParticle*) branchParticle->At(iGenB2);
      }
      else if ( deltaR(genParticle->Eta, vRecoB3.Eta(), genParticle->Phi, vRecoB3.Phi()) < MAX_MATCH_DIST ) {
	iGenB3 = iParticle;
	genB3 = (GenParticle*) branchParticle->At(iGenB3);
      }
      else if ( deltaR(genParticle->Eta, vRecoB4.Eta(), genParticle->Phi, vRecoB4.Phi()) < MAX_MATCH_DIST ) {
	iGenB4 = iParticle;
	genB4 = (GenParticle*) branchParticle->At(iGenB4);
      }
    }
    /*    cout << "Gen Particles " << endl;
    if (genB1) cout << "1 " << genB1->PT << " " << genB1->Eta << endl;
    if (genB2) cout << "2 " << genB2->PT << " " << genB2->Eta << endl;
    if (genB3) cout << "3 " << genB3->PT << " " << genB3->Eta << endl;
    if (genB4) cout << "4 " << genB4->PT << " " << genB4->Eta << endl;
    */

    for (Int_t iParticle=0; iParticle<branchParticle->GetEntries(); iParticle++) { // generator particle loop
      genParticle = (GenParticle*) branchParticle->At(iParticle);
      
      if (fabs(genParticle->PID) != H_ID_CODE) continue;

      if (iH1==-1) {
	iH1 = iParticle;
	genH1 = (GenParticle*) branchParticle->At(iH1);
      }
      else if (iH2==-1) {
	iH2 = iParticle;
	genH2 = (GenParticle*) branchParticle->At(iH2);
      }
    }

    /*    if (genH1) cout << "H1 " << genH1->PT << " " << genH1->Eta << " " << genH1->Phi <<  endl;
    if (genH2) cout << "H2 " << genH2->PT << " " << genH2->Eta << " " << genH2->Phi << endl;
    cout << endl;*/

    LorentzVector vGenB1(0,0,0,0);
    if (genB1) {
      vGenB1.SetPt(genB1->PT);
      vGenB1.SetEta(genB1->Eta);
      vGenB1.SetPhi(genB1->Phi);
      vGenB1.SetM(genB1->Mass);
      sGenB1 = &vGenB1;
    }
    else sGenB1 = &nothing;

    LorentzVector vGenB2(0,0,0,0);
    if (genB2) {
      vGenB2.SetPt(genB2->PT);
      vGenB2.SetEta(genB2->Eta);
      vGenB2.SetPhi(genB2->Phi);
      vGenB2.SetM(genB2->Mass);
      sGenB2 = &vGenB2;
    }
    else sGenB2 = &nothing;

    LorentzVector vGenB3(0,0,0,0);
    if (genB3) {
      vGenB3.SetPt(genB3->PT);
      vGenB3.SetEta(genB3->Eta);
      vGenB3.SetPhi(genB3->Phi);
      vGenB3.SetM(genB3->Mass);
      sGenB3 = &vGenB3;
    }
    else sGenB3 = &nothing;

    LorentzVector vGenB4(0,0,0,0);
    if (genB4) {
      vGenB4.SetPt(genB4->PT);
      vGenB4.SetEta(genB4->Eta);
      vGenB4.SetPhi(genB4->Phi);
      vGenB4.SetM(genB4->Mass);
      sGenB4 = &vGenB4;
    }
    else sGenB4 = &nothing;

    LorentzVector vGenH1(0,0,0,0);
    if (genH1) {
      vGenH1.SetPt(genH1->PT);
      vGenH1.SetEta(genH1->Eta);
      vGenH1.SetPhi(genH1->Phi);
      vGenH1.SetM(genH1->Mass);
      sGenH1 = &vGenH1;
    }
    else sGenH1 = &nothing;

    LorentzVector vGenH2(0,0,0,0);
    if (genH2) {
      vGenH2.SetPt(genH2->PT);
      vGenH2.SetEta(genH2->Eta);
      vGenH2.SetPhi(genH2->Phi);
      vGenH2.SetM(genH2->Mass);
      sGenH2 = &vGenH2;
    }
    else sGenH2 = &nothing;

    LorentzVector vTestBB1=vGenB1+vGenB2;
    LorentzVector vTestBB2=vGenB3+vGenB4;

    //cout << "test1 " << vTestBB1.Pt() << " " << vTestBB1.Eta() << " " << vTestBB1.Phi() << endl;
    //cout << "test1 " << vTestBB2.Pt() << " " << vTestBB2.Eta() << " " << vTestBB2.Phi() << endl;
    cout << "---" << endl;
    if ( (deltaR(vTestBB1.Eta(), vGenH1.Eta(), vTestBB1.Phi(), vGenH1.Phi()) < MAX_MATCH_DIST) && (deltaR(vTestBB2.Eta(), vGenH2.Eta(), vTestBB2.Phi(), vGenH2.Phi()) < MAX_MATCH_DIST) ) {
      cout << "H matched " << endl;
      nM++;
      iHmatch1=1;
      iHmatch2=1;
      iHmatch3=2;
      iHmatch4=2;
    }
    else if ( (deltaR(vTestBB2.Eta(), vGenH1.Eta(), vTestBB2.Phi(), vGenH1.Phi()) < MAX_MATCH_DIST) && (deltaR(vTestBB1.Eta(), vGenH2.Eta(), vTestBB1.Phi(), vGenH2.Phi()) < MAX_MATCH_DIST) ) {
      cout << "H matched " << endl;
      nM++;
      iHmatch1=2;
      iHmatch2=2;
      iHmatch3=1;
      iHmatch4=1;
    }

    vTestBB1=vGenB1+vGenB3;
    vTestBB2=vGenB2+vGenB4;

    //cout << "test2 " << vTestBB1.Pt() << " " << vTestBB1.Eta() << " " << vTestBB1.Phi() << endl;
    //cout << "test2 " << vTestBB2.Pt() << " " << vTestBB2.Eta() << " " << vTestBB2.Phi() << endl;

    if ( (deltaR(vTestBB1.Eta(), vGenH1.Eta(), vTestBB1.Phi(), vGenH1.Phi()) < MAX_MATCH_DIST) && (deltaR(vTestBB2.Eta(), vGenH2.Eta(), vTestBB2.Phi(), vGenH2.Phi()) < MAX_MATCH_DIST) ) {
      cout << "H matched " << endl;
      nM++;
      iHmatch1=1;
      iHmatch2=2;
      iHmatch3=1;
      iHmatch4=2;
    }
    else if ( (deltaR(vTestBB2.Eta(), vGenH1.Eta(), vTestBB2.Phi(), vGenH1.Phi()) < MAX_MATCH_DIST) && (deltaR(vTestBB1.Eta(), vGenH2.Eta(), vTestBB1.Phi(), vGenH2.Phi()) < MAX_MATCH_DIST) ) {
      cout << "H matched " << endl;
      nM++;
      iHmatch1=2;
      iHmatch2=1;
      iHmatch3=2;
      iHmatch4=1;
    }

    vTestBB1=vGenB1+vGenB4;
    vTestBB2=vGenB3+vGenB2;

    //cout << "test3 " << vTestBB1.Pt() << " " << vTestBB1.Eta() << " " << vTestBB1.Phi() << endl;
    //cout << "test3 " << vTestBB2.Pt() << " " << vTestBB2.Eta() << " " << vTestBB2.Phi() << endl;

    if ( (deltaR(vTestBB1.Eta(), vGenH1.Eta(), vTestBB1.Phi(), vGenH1.Phi()) < MAX_MATCH_DIST) && (deltaR(vTestBB2.Eta(), vGenH2.Eta(), vTestBB2.Phi(), vGenH2.Phi()) < MAX_MATCH_DIST) ) {
      cout << "H matched " << endl;
      nM++;
      iHmatch1=1;
      iHmatch2=2;
      iHmatch3=2;
      iHmatch4=1;
    }
    else if ( (deltaR(vTestBB2.Eta(), vGenH1.Eta(), vTestBB2.Phi(), vGenH1.Phi()) < MAX_MATCH_DIST) && (deltaR(vTestBB1.Eta(), vGenH2.Eta(), vTestBB1.Phi(), vGenH2.Phi()) < MAX_MATCH_DIST) ) {
      cout << "H matched " << endl;
      nM++;
      iHmatch1=2;
      iHmatch2=1;
      iHmatch3=1;
      iHmatch4=2;
    }

    if (iHmatch1==-1) {
      iHmatch1=0;
      iHmatch2=0;
      iHmatch3=0;
      iHmatch4=0;
    }
    //cout << iHmatch1 << " " << iHmatch2 << " " << iHmatch3 << " " << iHmatch4 << endl;
    //cout << endl;

    if (nM==0) iNot++;
    else if (nM==1) iMatched++;
    else if (nM>1) iTwo++;

    // match gen jets to reco jets
    for (Int_t iJet=0; iJet<branchGenJet->GetEntries(); iJet++) { // generator level jet loop
      genJet = (Jet*) branchGenJet->At(iJet);

      if ((jetB1) && (deltaR(genJet->Eta, jetB1->Eta, genJet->Phi, jetB1->Phi) < MAX_MATCH_DIST) ) {
	iGenJetB1=iJet;
	genJetB1 = (Jet*) branchGenJet->At(iGenJetB1);
      }
      else if ((jetB2) && (deltaR(genJet->Eta, jetB2->Eta, genJet->Phi, jetB2->Phi) < MAX_MATCH_DIST) ) {
	iGenJetB2=iJet;
	genJetB2 = (Jet*) branchGenJet->At(iGenJetB2);
      }
      else if ((jetB3) && (deltaR(genJet->Eta, jetB3->Eta, genJet->Phi, jetB3->Phi) < MAX_MATCH_DIST) ) {
	iGenJetB3=iJet;
	genJetB3 = (Jet*) branchGenJet->At(iGenJetB3);
      }
      else if ((jetB4) && (deltaR(genJet->Eta, jetB4->Eta, genJet->Phi, jetB4->Phi) < MAX_MATCH_DIST) ) {
	iGenJetB4=iJet;
	genJetB4 = (Jet*) branchGenJet->At(iGenJetB4);
      }
      else if ((jet1) && (deltaR(genJet->Eta, jet1->Eta, genJet->Phi, jet1->Phi) < MAX_MATCH_DIST) ) {
	iGenJetVBF1=iJet;
	genJetVBF1 = (Jet*) branchGenJet->At(iGenJetVBF1);
      }
      else if ((jet2) && (deltaR(genJet->Eta, jet2->Eta, genJet->Phi, jet2->Phi) < MAX_MATCH_DIST) ) {
	iGenJetVBF2=iJet;
	genJetVBF2 = (Jet*) branchGenJet->At(iGenJetVBF2);
      }
    }

    LorentzVector vGenJetB1(0,0,0,0);
    if (genJetB1) {
      vGenJetB1.SetPt(genJetB1->PT);
      vGenJetB1.SetEta(genJetB1->Eta);
      vGenJetB1.SetPhi(genJetB1->Phi);
      vGenJetB1.SetM(genJetB1->Mass);
      sGenJetB1 = &vGenJetB1;
    }
    else sGenJetB1 = &nothing;

    LorentzVector vGenJetB2(0,0,0,0);
    if (genJetB2) {
      vGenJetB2.SetPt(genJetB2->PT);
      vGenJetB2.SetEta(genJetB2->Eta);
      vGenJetB2.SetPhi(genJetB2->Phi);
      vGenJetB2.SetM(genJetB2->Mass);
      sGenJetB2 = &vGenJetB2;
    }
    else sGenJetB2 = &nothing;

    LorentzVector vGenJetB3(0,0,0,0);
    if (genJetB3) {
      vGenJetB3.SetPt(genJetB3->PT);
      vGenJetB3.SetEta(genJetB3->Eta);
      vGenJetB3.SetPhi(genJetB3->Phi);
      vGenJetB3.SetM(genJetB3->Mass);
      sGenJetB3 = &vGenJetB3;
    }
    else sGenJetB3 = &nothing;

    LorentzVector vGenJetB4(0,0,0,0);
    if (genJetB4) {
      vGenJetB4.SetPt(genJetB4->PT);
      vGenJetB4.SetEta(genJetB4->Eta);
      vGenJetB4.SetPhi(genJetB4->Phi);
      vGenJetB4.SetM(genJetB4->Mass);
      sGenJetB4 = &vGenJetB4;
    }
    else sGenJetB4 = &nothing;

    LorentzVector vGenJetVBF1(0,0,0,0);
    if (genJetVBF1) {
      vGenJetVBF1.SetPt(genJetVBF1->PT);
      vGenJetVBF1.SetEta(genJetVBF1->Eta);
      vGenJetVBF1.SetPhi(genJetVBF1->Phi);
      vGenJetVBF1.SetM(genJetVBF1->Mass);
      sGenJetVBF1 = &vGenJetVBF1;
    }
    else sGenJetVBF1 = &nothing;

    LorentzVector vGenJetVBF2(0,0,0,0);
    if (genJetVBF2) {
      vGenJetVBF2.SetPt(genJetVBF2->PT);
      vGenJetVBF2.SetEta(genJetVBF2->Eta);
      vGenJetVBF2.SetPhi(genJetVBF2->Phi);
      vGenJetVBF2.SetM(genJetVBF2->Mass);
      sGenJetVBF2 = &vGenJetVBF2;
    }
    else sGenJetVBF2 = &nothing;
    /*
    cout << "Gen Jets " << endl;
    if (genJetB1) cout << "1 " << genJetB1->PT << " " << genJetB1->Eta << endl;
    if (genJetB2) cout << "2 " << genJetB2->PT << " " << genJetB2->Eta << endl;
    if (genJetB3) cout << "3 " << genJetB3->PT << " " << genJetB3->Eta << endl;
    if (genJetB4) cout << "4 " << genJetB4->PT << " " << genJetB4->Eta << endl;
    if (genJetVBF1) cout << "V1 " << genJetVBF1->PT << " " << genJetVBF1->Eta << endl;
    if (genJetVBF2) cout << "V2 " << genJetVBF2->PT << " " << genJetVBF2->Eta << endl;
    cout << endl;
    */
    outTree->Fill();

  } // end event loop

  outFile->Write();
  outFile->Save();

  cout << endl;
  cout << "matched  : " << iMatched << endl;
  cout << "not      : " << iNot << endl;
  cout << "too much : " << iTwo << endl;
  cout << "total    : " << iMatched+iNot+iTwo << endl;

}

Float_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 ) {

  const Float_t pi = 3.14159265358979;

  Float_t etaDiff = (eta1-eta2);
  Float_t phiDiff = fabs(phi1-phi2);
  while ( phiDiff>pi ) phiDiff = fabs(phiDiff-2.0*pi);

  Float_t deltaRSquared = etaDiff*etaDiff + phiDiff*phiDiff;

  return TMath::Sqrt(deltaRSquared);

}

Int_t puJetID( Float_t eta, Float_t meanSqDeltaR, Float_t betastar) {
  
  Float_t MeanSqDeltaRMaxBarrel=0.07;
  Float_t BetaMinBarrel=0.87;
  Float_t MeanSqDeltaRMaxEndcap=0.07;
  Float_t BetaMinEndcap=0.85;

  //cout << eta << ", " << meanSqDeltaR << ", " << betastar << ": ";

  if (fabs(eta)<1.5) {
    if ((meanSqDeltaR<MeanSqDeltaRMaxBarrel)&&(betastar<BetaMinBarrel)) {
      //cout << "barrel 0" << endl;
      return 0;
    }
    else {
      //cout << "barrel 1" << endl;
      return 1;
    }
  }
  else if (fabs(eta)<4.0) {
    if ((meanSqDeltaR<MeanSqDeltaRMaxEndcap)&&(betastar<BetaMinEndcap)) {
      //cout << "endcap 0" << endl;
      return 0;
    }
    else {
      //cout << "endcap 1" << endl;
      return 1;
    }
  }
  //cout << "forward 1" << endl;
  return 1;

}
