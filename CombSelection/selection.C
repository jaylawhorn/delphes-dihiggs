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

//void selection(const TString inputfile="/afs/cern.ch/work/j/jlawhorn/public/new_hh/new_hh_0.root",
void selection(const TString inputfile="root://eoscms.cern.ch//store/group/upgrade/delphes/ProdJun14/tt-4p-1100-1700-v1510_14TEV/tt-4p-1100-1700-v1510_14TEV_151887068_PhaseII_Conf4_140PileUp_seed151887069_1of5.root",
	       const Float_t xsec=2.92,
	       const Float_t totalEvents=5000,
	       Int_t   sampleNo=100,
	       const TString outputfile="test.root") {

  cout << inputfile << " " << xsec << " " << totalEvents << " " << outputfile << endl;

  // declare constants
  const Double_t MUON_MASS = 0.105658369;
  const Double_t ELE_MASS  = 0.000511;

  const Int_t TAU_ID_CODE = 15;
  const Int_t B_ID_CODE = 5;
  const Int_t G_ID_CODE = 21;

  const Float_t MAX_MATCH_DIST = 0.4;

  // event categories
  enum { HH=0, H, TT, WJET, ZJET, EWK, ETC };

  // tau decay modes
  enum { hadron=1, electron, muon };

  // setup mt2 minimizer
  ROOT::Math::Minimizer *min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
  min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
  min->SetTolerance(10.0);
  min->SetPrintLevel(0);

  TVector2 tau1(0,0), tau2(0,0), mpt(0,0), ppmpt(0,0);
  TVector2 b1(0,0), b2(0,0);
  Float_t mTau1=0, mTau2=0;
  Float_t mB1=0, mB2=0;
  Double_t mt2=0;
  Double_t ppMt2=0;

  // read input input file
  TChain chain("Delphes");
  chain.Add(inputfile);
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  cout << numberOfEntries << endl; 
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
  TClonesArray *branchMET =treeReader->UseBranch("MissingET");
  TClonesArray *branchPuppiMET =treeReader->UseBranch("PuppiMissingET");

  TClonesArray *branchGenJet = treeReader->UseBranch("GenJet");
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchEvent = treeReader->UseBranch("Event");

  //set up loop variables
  GenParticle *genParticle=0;
  Jet *genJet=0;
  Jet *jet=0;
  Electron *ele=0;
  Muon *mu=0;
  Photon *gam=0;
  MissingET *missET=0;
  LHEFEvent *event=0;

  // set up storage variables
  // TAUS
  Jet *jetTau1=0, *jetTau2=0; Electron *eleTau=0; Muon *muTau=0;
  // B-JETS
  Jet *jetB1=0, *jetB2=0, *jetB3=0, *jetB4=0;
  // VBF-JETS
  Jet *jet_6j1=0, *jet_6j2=0, *jet_tt1=0, *jet_tt2=0;
  // PHOTONS
  Photon *gamma1=0, *gamma2=0;
  // GEN TAUS
  GenParticle *genTau1=0, *genTau2=0;
  // GEN B-QUARKS
  GenParticle *genB1=0, *genB2=0, *genB3=0, *genB4=0;
  // GEN PHOTONS
  GenParticle *genGam1=0, *genGam2=0;
  // HIGG(S)
  GenParticle *genH1=0, *genH2=0;
  // GEN B-JETS
  Jet *genJetB1=0, *genJetB2=0, *genJetB3=0, *genJetB4=0;
  // GEN TAU-JETS
  Jet *genJetTau1=0, *genJetTau2=0;
  // GEN VBF-JETS
  Jet *genJet_tt1=0, *genJet_tt2=0, *genJet_6j1=0, *genJet_6j2=0;

  Int_t iT1=-1, iT2=-1;
  Int_t iB1=-1, iB2=-1, iB3=-1, iB4=-1;
  Int_t iG1=-1, iG2=-1;
  Int_t iJ1=-1, iJ2=-1, iJ3=-1, iJ4=-1;
  Int_t iH1=-1, iH2=-1;

  Int_t iGenTau1=-1, iGenTau2=-1;
  Int_t iGenB1=-1,   iGenB2=-1, iGenB3=-1, iGenB4=-1;
  Int_t iGenGam1=-1, iGenGam2=-1;

  Int_t iGenJetTau1=-1, iGenJetTau2=-1;
  Int_t iGenJetB1=-1,   iGenJetB2=-1,   iGenJetB3=-1,   iGenJetB4=-1;
  Int_t iGenJet_tt1=-1, iGenJet_tt2=-1,  iGenJet_6j1=-1, iGenJet_6j2=-1;

  Int_t iHmatch1=-1, iHmatch2=-1, iHmatch3=-1, iHmatch4=-1;

  Int_t eventType;
  Int_t genInfo;
  Float_t eventWeight;

  Int_t isBBTT;
  Int_t isBBGG;
  Int_t isBBBB;
  Int_t isVBFTT;
  Int_t isVBF4B;

  Float_t met, metPhi;
  Float_t ppMet, ppMetPhi;

  Int_t nCentral=0, nBtag=0;
  Int_t centB=0;

  Int_t tauCat1=0, tauCat2=0;
  Int_t bTag1=0, bTag2=0, bTag3=0, bTag4=0;
  Int_t jbTag_tt1=0, jbTag_tt2=0, jbTag_6j1=0, jbTag_6j2=0;
  
  Float_t ptTau1, ptTau2, ptB1, ptB2, ptB3, ptB4, ptG1, ptG2, ptJet_tt1, ptJet_tt2, ptJet_6j1, ptJet_6j2;
  Float_t etaTau1, etaTau2, etaB1, etaB2, etaB3, etaB4, etaG1, etaG2, etaJet_tt1, etaJet_tt2, etaJet_6j1, etaJet_6j2;
  Float_t phiTau1, phiTau2, phiB1, phiB2, phiB3, phiB4, phiG1, phiG2, phiJet_tt1, phiJet_tt2, phiJet_6j1, phiJet_6j2;
  Float_t mB3, mB4, eG1, eG2, mJet_tt1, mJet_tt2, mJet_6j1, mJet_6j2;

  Float_t mTT=0, mBB=0, mGG=0;
  Float_t mHH=0, ptHH=0;
  Float_t mJJ_tt=0, dEta_tt=0;
  Float_t mJJ_6j=0, dEta_6j=0;

  LorentzVector *sRecoTau1=0,   *sRecoTau2=0;
  LorentzVector *sGenJetTau1=0, *sGenJetTau2=0;
  LorentzVector *sGenTau1=0,    *sGenTau2=0;

  LorentzVector *sRecoGam1=0,   *sRecoGam2=0;
  LorentzVector *sGenGam1=0,    *sGenGam2=0;

  LorentzVector *sRecoB1=0,   *sRecoB2=0,   *sRecoB3=0,   *sRecoB4=0;
  LorentzVector *sGenJetB1=0, *sGenJetB2=0, *sGenJetB3=0, *sGenJetB4=0;
  LorentzVector *sGenB1=0,    *sGenB2=0,    *sGenB3=0,    *sGenB4=0;

  LorentzVector *sRecoJet_tt1=0, *sRecoJet_tt2=0, *sRecoJet_6j1=0, *sRecoJet_6j2=0;
  LorentzVector *sGenJet_tt1=0,  *sGenJet_tt2=0,  *sGenJet_6j1=0,  *sGenJet_6j2=0;

  LorentzVector *sGenH1=0, *sGenH2=0;

  TFile *outFile = new TFile(outputfile, "RECREATE");

  // tree to hold information about selected events
  TTree *outTree = new TTree("Events", "Events");

  outTree->Branch("eventWeight",    &eventWeight,    "eventWeight/f");  // event weight from cross-section and Event->Weight
  outTree->Branch("sampleNo",       &sampleNo,       "sampleNo/i");     // sample number (see config card for details)
  outTree->Branch("genInfo",        &genInfo,        "genInfo/i");      // generator level info (see below)
  outTree->Branch("isBBTT",         &isBBTT,         "isBBTT/i");       // 1 if final state matches bbtt
  outTree->Branch("isBBGG",         &isBBGG,         "isBBGG/i");       // 1 if final state matches bbgg
  outTree->Branch("isBBBB",         &isBBBB,         "isBBBB/i");       // 1 if final state matches bbb
  outTree->Branch("isVBFTT",        &isVBFTT,        "isVBFTT/i");      // 1 if final state matches vbf htt
  outTree->Branch("isVBF4B",        &isVBF4B,        "isVBF4B/i");      // 1 if final state matches vbf hhbbbb
  outTree->Branch("eventType",      &eventType,      "eventType/i");    // see line 47

  outTree->Branch("met",            &met,            "met/f");          // missing transverse energy
  outTree->Branch("metPhi",         &metPhi,         "metPhi/f");       // missing transverse energy phi
  outTree->Branch("ppMet",          &ppMet,          "ppMet/f");        // PUPPI missing transverse energy
  outTree->Branch("ppMetPhi",       &ppMetPhi,       "ppMetPhi/f");     // PUPPI missing transverse energy phi

  outTree->Branch("ptTau1",         &ptTau1,         "ptTau1/f");       // pt(Tau1)
  outTree->Branch("etaTau1",        &etaTau1,        "etaTau1/f");      // eta(Tau1)
  outTree->Branch("phiTau1",        &phiTau1,        "phiTau1/f");      // phi(Tau1)
  outTree->Branch("mTau1",          &mTau1,          "mTau1/f");        // m(Tau1)
  outTree->Branch("tauCat1",        &tauCat1,        "tauCat1/i");      // leading tau final state - jet, muon, electron

  outTree->Branch("ptTau2",         &ptTau2,         "ptTau2/f");       // pt(Tau2)
  outTree->Branch("etaTau2",        &etaTau2,        "etaTau2/f");      // eta(Tau2)
  outTree->Branch("phiTau2",        &phiTau2,        "phiTau2/f");      // phi(Tau2)
  outTree->Branch("mTau2",          &mTau2,          "mTau2/f");        // m(Tau2)
  outTree->Branch("tauCat2",        &tauCat2,        "tauCat2/i");      // second tau final state - jet, muon, electron

  outTree->Branch("ptG1",           &ptG1,           "ptG1/f");         // pt(Gam1)
  outTree->Branch("etaG1",          &etaG1,          "etaG1/f");        // eta(Gam1)
  outTree->Branch("phiG1",          &phiG1,          "phiG1/f");        // phi(Gam1)
  outTree->Branch("eG1",            &eG1,            "eG1/f");          // e(Gam1)

  outTree->Branch("ptG2",           &ptG2,           "ptG2/f");         // pt(Gam2)
  outTree->Branch("etaG2",          &etaG2,          "etaG2/f");        // eta(Gam2)
  outTree->Branch("phiG2",          &phiG2,          "phiG2/f");        // phi(Gam2)
  outTree->Branch("eG2",            &eG2,            "eG2/f");          // e(Gam2)

  outTree->Branch("ptB1",           &ptB1,           "ptB1/f");         // pt(B1)
  outTree->Branch("etaB1",          &etaB1,          "etaB1/f");        // eta(B1)
  outTree->Branch("phiB1",          &phiB1,          "phiB1/f");        // phi(B1)
  outTree->Branch("mB1",            &mB1,            "mB1/f");          // m(B1)
  outTree->Branch("bTag1",          &bTag1,          "bTag1/i");        // leading b-jet tag from delphes
  outTree->Branch("iHmatch1",       &iHmatch1,       "iHmatch1/i");     // if 4 b's, which of two higgs matched

  outTree->Branch("ptB2",           &ptB2,           "ptB2/f");         // pt(B2)
  outTree->Branch("etaB2",          &etaB2,          "etaB2/f");        // eta(B2)
  outTree->Branch("phiB2",          &phiB2,          "phiB2/f");        // phi(B2)
  outTree->Branch("mB2",            &mB2,            "mB2/f");          // m(B2)
  outTree->Branch("bTag2",          &bTag2,          "bTag2/i");        // second b-jet tag from delphes
  outTree->Branch("iHmatch2",       &iHmatch2,       "iHmatch2/i");     // if 4 b's, which of two higgs matched

  outTree->Branch("ptB3",           &ptB3,           "ptB3/f");         // pt(B3)
  outTree->Branch("etaB3",          &etaB3,          "etaB3/f");        // eta(B3)
  outTree->Branch("phiB3",          &phiB3,          "phiB3/f");        // phi(B3)
  outTree->Branch("mB3",            &mB3,            "mB3/f");          // m(B3)
  outTree->Branch("bTag3",          &bTag3,          "bTag3/i");        // third b-jet tag from delphes
  outTree->Branch("iHmatch3",       &iHmatch3,       "iHmatch3/i");     // if 4 b's, which of two higgs matched

  outTree->Branch("ptB4",           &ptB4,           "ptB4/f");         // pt(B4)
  outTree->Branch("etaB4",          &etaB4,          "etaB4/f");        // eta(B4)
  outTree->Branch("phiB4",          &phiB4,          "phiB4/f");        // phi(B4)
  outTree->Branch("mB4",            &mB4,            "mB4/f");          // m(B4)
  outTree->Branch("bTag4",          &bTag4,          "bTag4/i");        // fourth b-jet tag from delphes
  outTree->Branch("iHmatch4",       &iHmatch4,       "iHmatch4/i");     // if 4 b's, which of two higgs matched

  outTree->Branch("ptJet_tt1",      &ptJet_tt1,      "ptJet_tt1/f");    // pt(Jet1)
  outTree->Branch("etaJet_tt1",     &etaJet_tt1,     "etaJet_tt1/f");   // eta(Jet1)
  outTree->Branch("phiJet_tt1",     &phiJet_tt1,     "phiJet_tt1/f");   // phi(Jet1)
  outTree->Branch("mJet_tt1",       &mJet_tt1,       "mJet_tt1/f");     // m(Jet1)
  outTree->Branch("jbTag_tt1",      &jbTag_tt1,      "jbTag_tt1/i");    // leading VBF-jet b tag from delphes

  outTree->Branch("ptJet_tt2",      &ptJet_tt2,      "ptJet_tt2/f");    // pt(Jet_tt2)
  outTree->Branch("etaJet_tt2",     &etaJet_tt2,     "etaJet_tt2/f");   // eta(Jet_tt2)
  outTree->Branch("phiJet_tt2",     &phiJet_tt2,     "phiJet_tt2/f");   // phi(Jet_tt2)
  outTree->Branch("mJet_tt2",       &mJet_tt2,       "mJet_tt2/f");     // m(Jet_tt2)
  outTree->Branch("jbTag_tt2",      &jbTag_tt2,      "jbTag_tt2/i");    // second VBF-jet b tag from delphes

  outTree->Branch("ptJet_6j1",      &ptJet_6j1,      "ptJet_6j1/f");    // pt(Jet1)
  outTree->Branch("etaJet_6j1",     &etaJet_6j1,     "etaJet_6j1/f");   // eta(Jet1)
  outTree->Branch("phiJet_6j1",     &phiJet_6j1,     "phiJet_6j1/f");   // phi(Jet1)
  outTree->Branch("mJet_6j1",       &mJet_6j1,       "mJet_6j1/f");     // m(Jet1)
  outTree->Branch("jbTag_6j1",      &jbTag_6j1,      "jbTag_6j1/i");    // leading VBF-jet b tag from delphes

  outTree->Branch("ptJet_6j2",      &ptJet_6j2,      "ptJet_6j2/f");    // pt(Jet2)
  outTree->Branch("etaJet_6j2",     &etaJet_6j2,     "etaJet_6j2/f");   // eta(Jet2)
  outTree->Branch("phiJet_6j2",     &phiJet_6j2,     "phiJet_6j2/f");   // phi(Jet2)
  outTree->Branch("mJet_6j2",       &mJet_6j2,       "mJet_6j2/f");     // m(Jet2)
  outTree->Branch("jbTag_6j2",      &jbTag_6j2,      "jbTag_6j2/i");    // second VBF-jet b tag from delphes

  outTree->Branch("mTT",            &mTT,            "mTT/f");          // mass(tautau)
  outTree->Branch("mGG",            &mGG,            "mGG/f");          // mass(gamgam)
  outTree->Branch("mBB",            &mBB,            "mBB/f");          // mass(bb)

  outTree->Branch("mt2",            &mt2,            "mt2/D");          // "stransverse mass" (HH)
  outTree->Branch("ppMt2",          &ppMt2,          "ppMt2/D");        // PUPPI "stransverse mass" (HH)
  outTree->Branch("mHH",            &mHH,            "mHH/f");          // mass(HH) 
  outTree->Branch("ptHH",           &ptHH,           "ptHH/f");         // pt(HH)

  outTree->Branch("nBtag",          &nBtag,          "nBtag/i");        // number of b-tagged jets (VBF)   
  outTree->Branch("nCentral",       &nCentral,       "nCentral/i");     // number of central jets (VBF)
  outTree->Branch("centB",          &centB,          "centB/i");        // how many b's between VBF jets?
  outTree->Branch("mJJ_tt",         &mJJ_tt,         "mJJ_tt/f");       // mass(jetjet) (VBF)
  outTree->Branch("dEta_tt",        &dEta_tt,        "dEta_tt/f");      // delta Eta (VBF)
  outTree->Branch("mJJ_6j",         &mJJ_6j,         "mJJ_6j/f");       // mass(jetjet) (VBF)
  outTree->Branch("dEta_6j",        &dEta_6j,        "dEta_6j/f");      // delta Eta (VBF)

  outTree->Branch("sGenTau1",       "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sGenTau1);      // 4-vector for generator leading tau
  outTree->Branch("sGenJetTau1",    "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sGenJetTau1);   // 4-vector for generator leading tau
  outTree->Branch("sRecoTau1",      "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sRecoTau1);     // 4-vector for reconstructed leading tau

  outTree->Branch("sGenTau2",       "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sGenTau2);      // 4-vector for generator second tau
  outTree->Branch("sGenJetTau2",    "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sGenJetTau2);   // 4-vector for generator second tau
  outTree->Branch("sRecoTau2",      "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sRecoTau2);     // 4-vector for reconstructed second tau

  outTree->Branch("sGenGam1",       "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sGenGam1);      // 4-vector for generator leading gamma
  outTree->Branch("sRecoGam1",      "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sRecoGam1);     // 4-vector for reconstructed leading gamma

  outTree->Branch("sGenGam2",       "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sGenGam2);      // 4-vector for generator second gamma
  outTree->Branch("sRecoGam2",      "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sRecoGam2);     // 4-vector for reconstructed second gamma

  outTree->Branch("sGenB1",         "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sGenB1);        // 4-vector for generator leading b-quark
  outTree->Branch("sGenJetB1",      "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sGenJetB1);     // 4-vector for generator leading b-jet
  outTree->Branch("sRecoB1",        "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sRecoB1);       // 4-vector for reconstructed leading b-jet

  outTree->Branch("sGenB2",         "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sGenB2);        // 4-vector for generator b-quark
  outTree->Branch("sGenJetB2",      "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sGenJetB2);     // 4-vector for generator b-jet
  outTree->Branch("sRecoB2",        "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sRecoB2);       // 4-vector for reconstructed b-jet

  outTree->Branch("sGenB3",         "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sGenB3);        // 4-vector for generator b-quark
  outTree->Branch("sGenJetB3",      "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sGenJetB3);     // 4-vector for generator b-jet
  outTree->Branch("sRecoB3",        "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sRecoB3);       // 4-vector for reconstructed b-jet

  outTree->Branch("sGenB4",         "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sGenB4);        // 4-vector for generator b-quark
  outTree->Branch("sGenJetB4",      "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sGenJetB4);     // 4-vector for generator b-jet
  outTree->Branch("sRecoB4",        "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sRecoB4);       // 4-vector for reconstructed b-jet

  outTree->Branch("sRecoJet_tt1",   "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sRecoJet_tt1);  // 4-vector for reconstructed leading VBF-jet
  outTree->Branch("sGenJet_tt1",    "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sGenJet_tt1);   // 4-vector for generator leading VBF-jet

  outTree->Branch("sRecoJet_tt2",   "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sRecoJet_tt2);  // 4-vector for reconstructed VBF-jet
  outTree->Branch("sGenJet_tt2",    "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sGenJet_tt2);   // 4-vector for generator VBF-jet

  outTree->Branch("sRecoJet_6j1",   "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sRecoJet_6j1);  // 4-vector for reconstructed leading VBF-jet
  outTree->Branch("sGenJet_6j1",    "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sGenJet_6j1);   // 4-vector for generator leading VBF-jet

  outTree->Branch("sRecoJet_6j2",   "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sRecoJet_6j2);  // 4-vector for reconstructed VBF-jet
  outTree->Branch("sGenJet_6j2",    "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sGenJet_6j2);   // 4-vector for generator VBF-jet

  outTree->Branch("sGenH1",         "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sGenH1);        // generator Higgs
  outTree->Branch("sGenH2",         "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sGenH2);        // generator Higgs

  // define placeholder vector for things that don't exist
  LorentzVector nothing(-999,-999,0,-999);

  for (Int_t iEntry=0; iEntry<numberOfEntries; iEntry++) { // entry loop
    treeReader->ReadEntry(iEntry);

    // ********************
    // RESET
    // ********************

    iB1=-1; iB2=-1; iB3=-1; iB4=-1; iG1=-1; iG2=-1; 
    iT1=-1; iT2=-1; iJ1=-1; iJ2=-1; iJ3=-1; iJ4=-1;
    iH1=-1; iH2=-1;

    iGenTau1=-1;    iGenTau2=-1;
    iGenB1=-1;      iGenB2=-1;
    iGenB3=-1;      iGenB4=-1;
    iGenGam1=-1;    iGenGam2=-1;

    iGenJetTau1=-1; iGenJetTau2=-1;
    iGenJetB1=-1;   iGenJetB2=-1;
    iGenJetB3=-1;   iGenJetB4=-1;
    iGenJet_tt1=-1; iGenJet_tt2=-1;
    iGenJet_6j1=-1; iGenJet_6j2=-1;

    met=0; metPhi=0; ppMet=0; ppMetPhi=0;

    tauCat1=0; tauCat2=0; 
    bTag1=0; bTag2=0; bTag3=0; bTag4=0;
    jbTag_tt1=0; jbTag_tt2=0;
    jbTag_6j1=0; jbTag_6j2=0;
    
    mTT=-999; mBB=-999; mGG=-999;
    mHH=-999; ptHH=-999;
    mJJ_tt=-999; dEta_tt=-999; 
    mJJ_6j=-999; dEta_6j=-999; 
    nBtag=0; nCentral=0; centB=0;

    ptTau1=-999; etaTau1=-999; phiTau1=-999; mTau1=-999;
    ptTau2=-999; etaTau2=-999; phiTau2=-999; mTau2=-999;
    ptG1=-999; etaG1=-999; phiG1=-999; eG1=-999;
    ptG2=-999; etaG2=-999; phiG2=-999; eG2=-999;
    ptB1=-999; etaB1=-999; phiB1=-999; mB1=-999;
    ptB2=-999; etaB2=-999; phiB2=-999; mB2=-999;
    ptB3=-999; etaB3=-999; phiB3=-999; mB3=-999;
    ptB4=-999; etaB4=-999; phiB4=-999; mB4=-999;

    eventType=-1;

    iHmatch1=0; iHmatch2=0; iHmatch3=0; iHmatch4=0;

    isBBTT=0; isBBGG=0; isBBBB=0;
    isVBFTT=0; isVBF4B=0;

    jetTau1=0; jetTau2=0; eleTau=0;  muTau=0; 
    jetB1=0;   jetB2=0;   jetB3=0;   jetB4=0;
    jet_tt1=0; jet_tt2=0; jet_6j1=0; jet_6j1=0;
    gamma1=0;  gamma2=0;

    genTau1=0; genTau2=0; genGam1=0; genGam2=0;
    genB1=0;   genB2=0;   genB3=0;   genB4=0;

    genJetTau1=0; genJetTau2=0; 
    genJet_tt1=0; genJet_tt2=0; genJet_6j1=0; genJet_6j2=0; 
    genJetB1=0;   genJetB2=0;   genJetB3=0;   genJetB4=0; 

    sRecoB1=0;   sRecoB2=0;   sRecoB3=0;     sRecoB4=0;
    sRecoTau1=0; sRecoTau2=0; 
    sRecoJet_tt1=0; sRecoJet_tt2=0; sRecoJet_6j1=0; sRecoJet_6j2=0;
    sRecoGam1=0; sRecoGam2=0;
 
    sGenJetB1=0;   sGenJetB2=0;   sGenJetB3=0;   sGenJetB4=0;
    sGenJetTau1=0; sGenJetTau2=0; 
    sGenJet_tt1=0; sGenJet_tt2=0; sGenJet_6j1=0; sGenJet_6j2=0;

    sGenB1=0;   sGenB2=0;   sGenB3=0;   sGenB4=0;
    sGenTau1=0; sGenTau2=0; sGenGam1=0; sGenGam2=0;
    sGenH1=0;   sGenH2=0;

    // ********************
    // EVENT WEIGHT
    // ********************

    eventWeight = 1;
    if (branchEvent) {
      event = (LHEFEvent*) branchEvent->At(0);
      eventWeight*=event->Weight;
    }
    eventWeight *= xsec;
    eventWeight /= totalEvents;

    // ********************
    // TAU SELECTION
    // ********************
    
    for (Int_t iJet=0; iJet<branchJet->GetEntries(); iJet++) { // reconstructed jet loop
      jet = (Jet*) branchJet->At(iJet);
      
      if (fabs(jet->Eta)>4.0) continue;
      if (jet->PT<30) continue;
      if (jet->TauTag==0) continue;
      
      if (puJetID(jet->Eta, jet->MeanSqDeltaR, jet->BetaStar)==1) continue;

      if (iT1==-1) { 
	iT1=iJet; 
	jetTau1 = (Jet*) branchJet->At(iT1); 
	tauCat1=hadron;
      }
      else if (jet->PT > jetTau1->PT) { 
	iT2=iT1; 
	jetTau2 = (Jet*) branchJet->At(iT2); 
	tauCat2=hadron; 
	iT1=iJet; 
	jetTau1 = (Jet*) branchJet->At(iT1); 
	tauCat1=hadron;
      }
      else if (iT2==-1) { 
	iT2=iJet; 
	jetTau2 = (Jet*) branchJet->At(iT2); 
	tauCat2=hadron;
      }
      else if (jet->PT > jetTau2->PT) { 
	iT2=iJet; 
	jetTau2 = (Jet*) branchJet->At(iT2); 
	tauCat2=hadron;
      }
    } // end reco jet loop

    if ((iT1==-1) || (iT2==-1)) {
      // get muonic taus
      for (Int_t iMuon=0; iMuon<branchMuon->GetEntries(); iMuon++) { // reco muon loop
	mu = (Muon*) branchMuon->At(iMuon);

	if (fabs(mu->Eta)>4.0) continue;
	if (mu->PT<30) continue;

	if (iT1==-1 && tauCat2!=muon) { 
	  iT1=iMuon; 
	  muTau = (Muon*) branchMuon->At(iT1); 
	  tauCat1=muon; 
	}
	else if (iT2==-1 && tauCat1!=muon) { 
	  iT2=iMuon; 
	  muTau = (Muon*) branchMuon->At(iT2); 
	  tauCat2=muon; 
	}
	else if (muTau) {
	  if ( mu->PT > muTau->PT ) { 
	    if (tauCat1==muon) {
	      iT1=iMuon;
	      muTau = (Muon*) branchMuon->At(iT1); 
	    }
	    else if (tauCat2==muon) {
	      iT2=iMuon;
	      muTau = (Muon*) branchMuon->At(iT2); 
	    }
	  }
	}
      }

      // get electronic taus
      for (Int_t iEle=0; iEle<branchElectron->GetEntries(); iEle++) { // reco ele loop
	ele = (Electron*) branchElectron->At(iEle);

	if (fabs(ele->Eta)>4.0) continue;
	if (ele->PT<30) continue;

	if (iT1==-1 && tauCat1!=electron) { 
	  iT1=iEle; 
	  eleTau = (Electron*) branchElectron->At(iT1); 
	  tauCat1=electron; 
	}
	else if (iT2==-1 && tauCat2!=electron) { 
	  iT2=iEle; 
	  eleTau = (Electron*) branchElectron->At(iT2); 
	  tauCat2=electron;
	}
	else if (eleTau) {
	  if ( ele->PT > eleTau->PT ) { 
	    if ( tauCat1==electron) {
	      iT1=iEle;
	      eleTau = (Electron*) branchElectron->At(iT1); 
	    }
	    else if ( tauCat2==electron) {
	      iT2=iEle;
	      eleTau = (Electron*) branchElectron->At(iT2); 
	    }
	  }
	}
      }
    }

    // ********************
    // GAMMA SELECTION
    // ********************

    for (Int_t iP=0; iP<branchPhoton->GetEntries(); iP++) { // reco photon loop
      gam = (Photon*) branchPhoton->At(iP);

      if (fabs(gam->Eta)>4.0) continue;
      if (gam->PT<30) continue;

      if (iG1==-1) {
	iG1=iP;
	gamma1 = (Photon*) branchPhoton->At(iG1);
      }
      else if (gam->PT>gamma1->PT) {
	iG2=iG1;
	gamma2 = (Photon*) branchPhoton->At(iG2);
	iG1=iP;
	gamma1 = (Photon*) branchPhoton->At(iG1);
      }
      else if (iG2==-1) {
	iG2=iP;
	gamma2 = (Photon*) branchPhoton->At(iG2);
      }
      else if (gam->PT>gamma2->PT) {
	iG1=iP;
	gamma1 = (Photon*) branchPhoton->At(iG1);
      }
    } // end reco photon loop

    // ********************
    // B-JET SELECTION
    // ********************

    for (Int_t iJet=0; iJet<branchJet->GetEntries(); iJet++) { // reconstructed jet loop
      jet = (Jet*) branchJet->At(iJet);

      if (fabs(jet->Eta)>4.0) continue;
      if (jet->PT<30) continue;
      if (jet->BTag==0) continue;

      if (puJetID(jet->Eta, jet->MeanSqDeltaR, jet->BetaStar)==1) continue;

      if ((jetTau1)&&(deltaR(jet->Eta, jetTau1->Eta, jet->Phi, jetTau1->Phi) < MAX_MATCH_DIST)) continue;
      if ((jetTau2)&&(deltaR(jet->Eta, jetTau2->Eta, jet->Phi, jetTau2->Phi) < MAX_MATCH_DIST)) continue;

      if (iB1==-1) {
	iB1=iJet; 
	jetB1 = (Jet*) branchJet->At(iB1); 
	bTag1=jetB1->BTag;
      }
      else if (jet->PT > jetB1->PT) {
	if (jetB3) {
	  iB4=iB3;
	  jetB4 = (Jet*) branchJet->At(iB4);
	  bTag4=jetB4->BTag;
	}
	if (jetB2) {
	  iB3=iB2;
	  jetB3 = (Jet*) branchJet->At(iB3);
	  bTag3=jetB3->BTag;
	}
	iB2=iB1; 
	jetB2 = (Jet*) branchJet->At(iB2); 
	bTag2=jetB2->BTag;
	iB1=iJet;
	jetB1 = (Jet*) branchJet->At(iB1);
	bTag1=jetB1->BTag;
      }
      else if (iB2==-1) { 
	iB2=iJet; 
	jetB2 = (Jet*) branchJet->At(iB2); 
	bTag2=jetB2->BTag;
      }
      else if (jet->PT > jetB2->PT) { 
	if (jetB3) {
	  iB4=iB3;
	  jetB4 = (Jet*) branchJet->At(iB4);
	  bTag4=jetB4->BTag;
	}
	iB3=iB2;
	jetB3 = (Jet*) branchJet->At(iB3);
	bTag3=jetB3->BTag;
	iB2=iJet; 
	jetB2 = (Jet*) branchJet->At(iB2); 
	bTag2=jetB2->BTag;
      }
      else if (iB3==-1) {
	iB3=iJet;
	jetB3 = (Jet*) branchJet->At(iB3);
	bTag3=jetB3->BTag;
      }
      else if (jet->PT > jetB3->PT) {
	iB4=iB3;
	jetB4 = (Jet*) branchJet->At(iB4);
	bTag4=jetB4->BTag;
	iB3=iJet;
	jetB3 = (Jet*) branchJet->At(iB3);
	bTag3=jetB3->BTag;
      }
      else if (iB4==-1) {
	iB4=iJet;
	jetB4 = (Jet*) branchJet->At(iB4);
	bTag4=jetB4->BTag;
      }
      else if (jet->PT > jetB4->PT) {
	iB4=iJet;
	jetB4 = (Jet*) branchJet->At(iB4);
	bTag4=jetB4->BTag;
      }
    }

    // ********************
    // VBF-JET SELECTION
    // ********************

    for (Int_t iJet=0; iJet<branchJet->GetEntries(); iJet++) { // reconstructed jet loop
      jet = (Jet*) branchJet->At(iJet);

      if (jet->BTag>0) nBtag++;

      if (fabs(jet->Eta)>4.7) continue;
      if (jet->PT<30) continue;

      if (puJetID(jet->Eta, jet->MeanSqDeltaR, jet->BetaStar)==1) continue;

      if ((jetTau1)&&(deltaR(jet->Eta, jetTau1->Eta, jet->Phi, jetTau1->Phi) < MAX_MATCH_DIST)) continue;
      if ((jetTau2)&&(deltaR(jet->Eta, jetTau2->Eta, jet->Phi, jetTau2->Phi) < MAX_MATCH_DIST)) continue;
      
      if (iJ1==-1) {
        iJ1=iJet;
        jet_tt1 = (Jet*) branchJet->At(iJ1);
        jbTag_tt1=jet_tt1->BTag;
      }
      else if (jet->PT > jet_tt1->PT) {
        iJ2=iJ1;
        jet_tt2 = (Jet*) branchJet->At(iJ2);
        jbTag_tt2=jet_tt2->BTag;
        iJ1=iJet;
        jet_tt1 = (Jet*) branchJet->At(iJ1);
        jbTag_tt1=jet_tt1->BTag;
      }
      else if (iJ2==-1) {
        iJ2=iJet;
        jet_tt2 = (Jet*) branchJet->At(iJ2);
        jbTag_tt2=jet_tt2->BTag;
      }
      else if (jet->PT > jet_tt2->PT) {
        iJ2=iJet;
        jet_tt2 = (Jet*) branchJet->At(iJ2);
        jbTag_tt2=jet_tt2->BTag;
      }
    }

    for (Int_t iJet=0; iJet<branchJet->GetEntries(); iJet++) { // reconstructed jet loop
      jet = (Jet*) branchJet->At(iJet);

      if (fabs(jet->Eta)>4.7) continue;
      if (jet->PT<30) continue;

      if (puJetID(jet->Eta, jet->MeanSqDeltaR, jet->BetaStar)==1) continue;

      if ((jetB1)&&(deltaR(jet->Eta, jetB1->Eta, jet->Phi, jetB1->Phi) < MAX_MATCH_DIST)) continue;
      if ((jetB2)&&(deltaR(jet->Eta, jetB2->Eta, jet->Phi, jetB2->Phi) < MAX_MATCH_DIST)) continue;
      if ((jetB3)&&(deltaR(jet->Eta, jetB3->Eta, jet->Phi, jetB3->Phi) < MAX_MATCH_DIST)) continue;
      if ((jetB4)&&(deltaR(jet->Eta, jetB4->Eta, jet->Phi, jetB4->Phi) < MAX_MATCH_DIST)) continue;

      if (iJ3==-1) {
        iJ3=iJet;
        jet_6j1 = (Jet*) branchJet->At(iJ3);
        jbTag_6j1=jet_6j1->BTag;
      }
      else if (jet->PT > jet_6j1->PT) {
	iJ4=iJ3;
	jet_6j2 = (Jet*) branchJet->At(iJ4);
        jbTag_6j2=jet_6j2->BTag;
        iJ3=iJet;
        jet_6j1 = (Jet*) branchJet->At(iJ3);
        jbTag_6j1=jet_6j1->BTag;
      }
      else if (iJ4==-1) {
        iJ4=iJet;
        jet_6j2 = (Jet*) branchJet->At(iJ4);
        jbTag_6j2=jet_6j2->BTag;
      }
      else if (jet->PT > jet_6j2->PT) {
        iJ4=iJet;
        jet_6j2 = (Jet*) branchJet->At(iJ4);
        jbTag_6j2=jet_6j2->BTag;
      }
    }

    // ********************
    // STORE VARIABLES
    // ********************

    // fill 4-vector for leading tau
    LorentzVector vRecoTau1(0,0,0,0);
    if (jetTau1) {
      vRecoTau1.SetPt(jetTau1->PT);
      vRecoTau1.SetEta(jetTau1->Eta);
      vRecoTau1.SetPhi(jetTau1->Phi);
      vRecoTau1.SetM(jetTau1->Mass);
      sRecoTau1 = &vRecoTau1;
      ptTau1=jetTau1->PT;
      etaTau1=jetTau1->Eta;
      phiTau1=jetTau1->Phi;
      mTau1=jetTau1->Mass;
    }
    else if ((muTau)&&(tauCat1==muon)) {
      vRecoTau1.SetPt(muTau->PT);
      vRecoTau1.SetEta(muTau->Eta);
      vRecoTau1.SetPhi(muTau->Phi);
      vRecoTau1.SetM(MUON_MASS);
      sRecoTau1 = &vRecoTau1;
      ptTau1=muTau->PT;
      etaTau1=muTau->Eta;
      phiTau1=muTau->Phi;
      mTau1=MUON_MASS;
    }
    else if ((eleTau)&&(tauCat1==electron)) {
      vRecoTau1.SetPt(eleTau->PT);
      vRecoTau1.SetEta(eleTau->Eta);
      vRecoTau1.SetPhi(eleTau->Phi);
      vRecoTau1.SetM(ELE_MASS);
      sRecoTau1 = &vRecoTau1;
      ptTau1=eleTau->PT;
      etaTau1=eleTau->Eta;
      phiTau1=eleTau->Phi;
      mTau1=ELE_MASS;
    }
    else sRecoTau1 = &nothing;

    // fill 4-vector for second tau
    LorentzVector vRecoTau2(0,0,0,0);
    if (jetTau2) {
      vRecoTau2.SetPt(jetTau2->PT);
      vRecoTau2.SetEta(jetTau2->Eta);
      vRecoTau2.SetPhi(jetTau2->Phi);
      vRecoTau2.SetM(jetTau2->Mass);
      sRecoTau2 = &vRecoTau2;
      ptTau2=jetTau2->PT;
      etaTau2=jetTau2->Eta;
      phiTau2=jetTau2->Phi;
      mTau2=jetTau2->Mass;
    }
    else if ((muTau)&&(tauCat2==muon)) {
      vRecoTau2.SetPt(muTau->PT);
      vRecoTau2.SetEta(muTau->Eta);
      vRecoTau2.SetPhi(muTau->Phi);
      vRecoTau2.SetM(MUON_MASS);
      sRecoTau2 = &vRecoTau2;
      ptTau2=muTau->PT;
      etaTau2=muTau->Eta;
      phiTau2=muTau->Phi;
      mTau2=MUON_MASS;
    }
    else if ((eleTau)&&(tauCat2==electron)) {
      vRecoTau2.SetPt(eleTau->PT);
      vRecoTau2.SetEta(eleTau->Eta);
      vRecoTau2.SetPhi(eleTau->Phi);
      vRecoTau2.SetM(ELE_MASS);
      sRecoTau2 = &vRecoTau2;
      ptTau2=eleTau->PT;
      etaTau2=eleTau->Eta;
      phiTau2=eleTau->Phi;
      mTau2=ELE_MASS;
    }
    else sRecoTau2 = &nothing;

    // fill 4-vector for leading b-jet
    LorentzVector vRecoB1(0,0,0,0);
    if (jetB1) {
      vRecoB1.SetPt(jetB1->PT);
      vRecoB1.SetEta(jetB1->Eta);
      vRecoB1.SetPhi(jetB1->Phi);
      vRecoB1.SetM(jetB1->Mass);
      sRecoB1 = &vRecoB1;
      ptB1=jetB1->PT;
      etaB1=jetB1->Eta;
      phiB1=jetB1->Phi;
      mB1=jetB1->Mass;
    }
    else sRecoB1 = &nothing;

    // fill 4-vector for second b-jet
    LorentzVector vRecoB2(0,0,0,0);
    if (jetB2) {
      vRecoB2.SetPt(jetB2->PT);
      vRecoB2.SetEta(jetB2->Eta);
      vRecoB2.SetPhi(jetB2->Phi);
      vRecoB2.SetM(jetB2->Mass);
      sRecoB2 = &vRecoB2;
      ptB2=jetB2->PT;
      etaB2=jetB2->Eta;
      phiB2=jetB2->Phi;
      mB2=jetB2->Mass;
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
      ptB3=jetB3->PT;
      etaB3=jetB3->Eta;
      phiB3=jetB3->Phi;
      mB3=jetB3->Mass;
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
      ptB4=jetB4->PT;
      etaB4=jetB4->Eta;
      phiB4=jetB4->Phi;
      mB4=jetB4->Mass;
    }
    else sRecoB4 = &nothing;

    // fill 4-vector for leading photon
    LorentzVector vRecoGam1(0,0,0,0);
    if (gamma1) {
      vRecoGam1.SetPt(gamma1->PT);
      vRecoGam1.SetEta(gamma1->Eta);
      vRecoGam1.SetPhi(gamma1->Phi);
      vRecoGam1.SetM(0);
      sRecoGam1 = &vRecoGam1;
      ptG1=gamma1->PT;
      etaG1=gamma1->Eta;
      phiG1=gamma1->Phi;
      eG1=gamma1->E;
    }
    else sRecoGam1 = &nothing;

    // fill 4-vector for second photon
    LorentzVector vRecoGam2(0,0,0,0);
    if (gamma2) {
      vRecoGam2.SetPt(gamma2->PT);
      vRecoGam2.SetEta(gamma2->Eta);
      vRecoGam2.SetPhi(gamma2->Phi);
      vRecoGam2.SetM(0);
      sRecoGam2 = &vRecoGam2;
      ptG2=gamma2->PT;
      etaG2=gamma2->Eta;
      phiG2=gamma2->Phi;
      eG2=gamma2->E;
    }
    else sRecoGam2 = &nothing;

    // fill 4-vector for leading VBF jet
    LorentzVector vRecoJet_tt1(0,0,0,0);
    if (jet_tt1) {
      vRecoJet_tt1.SetPt(jet_tt1->PT);
      vRecoJet_tt1.SetEta(jet_tt1->Eta);
      vRecoJet_tt1.SetPhi(jet_tt1->Phi);
      vRecoJet_tt1.SetM(jet_tt1->Mass);
      sRecoJet_tt1 = &vRecoJet_tt1;
      ptJet_tt1=jet_tt1->PT;
      etaJet_tt1=jet_tt1->Eta;
      phiJet_tt1=jet_tt1->Phi;
      mJet_tt1=jet_tt1->Mass;
    }
    else sRecoJet_tt1 = &nothing;

    LorentzVector vRecoJet_tt2(0,0,0,0);
    if (jet_tt2) {
      vRecoJet_tt2.SetPt(jet_tt2->PT);
      vRecoJet_tt2.SetEta(jet_tt2->Eta);
      vRecoJet_tt2.SetPhi(jet_tt2->Phi);
      vRecoJet_tt2.SetM(jet_tt2->Mass);
      sRecoJet_tt2 = &vRecoJet_tt2;
      ptJet_tt2=jet_tt2->PT;
      etaJet_tt2=jet_tt2->Eta;
      phiJet_tt2=jet_tt2->Phi;
      mJet_tt2=jet_tt2->Mass;
    }
    else sRecoJet_tt2 = &nothing;

    // fill 4-vector for leading VBF jet
    LorentzVector vRecoJet_6j1(0,0,0,0);
    if (jet_6j1) {
      vRecoJet_6j1.SetPt(jet_6j1->PT);
      vRecoJet_6j1.SetEta(jet_6j1->Eta);
      vRecoJet_6j1.SetPhi(jet_6j1->Phi);
      vRecoJet_6j1.SetM(jet_6j1->Mass);
      sRecoJet_6j1 = &vRecoJet_6j1;
      ptJet_6j1=jet_6j1->PT;
      etaJet_6j1=jet_6j1->Eta;
      phiJet_6j1=jet_6j1->Phi;
      mJet_6j1=jet_6j1->Mass;
    }
    else sRecoJet_6j1 = &nothing;

    LorentzVector vRecoJet_6j2(0,0,0,0);
    if (jet_6j2) {
      vRecoJet_6j2.SetPt(jet_6j2->PT);
      vRecoJet_6j2.SetEta(jet_6j2->Eta);
      vRecoJet_6j2.SetPhi(jet_6j2->Phi);
      vRecoJet_6j2.SetM(jet_6j2->Mass);
      sRecoJet_6j2 = &vRecoJet_6j2;
      ptJet_6j2=jet_6j2->PT;
      etaJet_6j2=jet_6j2->Eta;
      phiJet_6j2=jet_6j2->Phi;
      mJet_6j2=jet_6j2->Mass;
    }
    else sRecoJet_6j2 = &nothing;
    
    // ********************
    // COMPUTE VARIABLES
    // ********************

    LorentzVector vTT;
    if (vRecoTau1.Pt()>0 && vRecoTau2.Pt()>0) {
      vTT = vRecoTau1+vRecoTau2;
      mTT=vTT.M();
    }
    else {
      mTT=-999;
      vTT=nothing;
    }
    
    LorentzVector vBB;
    if (vRecoB1.Pt()>0 && vRecoB2.Pt()>0) {
      vBB = vRecoB1+vRecoB2;
      mBB=vBB.M();
    }
    else {
      mBB=-999;
      vBB=nothing;
    }
    
    if (vRecoTau1.Pt()>0 && vRecoTau2.Pt()>0 && vRecoB1.Pt()>0 && vRecoB2.Pt()>0) {
      LorentzVector vHH = vTT+vBB;
      mHH=vHH.M();
      ptHH=vHH.Pt();
    }
    else {
      mHH=-999;
      ptHH=-999;
    }

    LorentzVector vGG;
    if (vRecoGam1.Pt()>0 && vRecoGam2.Pt()>0) {
      vGG = vRecoGam1+vRecoGam2;
      mGG=vGG.M();
    }
    else {
      mGG=-999;
      vGG=nothing;
    }

    if (vRecoGam1.Pt()>0 && vRecoGam2.Pt()>0 && vRecoB1.Pt()>0 && vRecoB2.Pt()>0) {
      LorentzVector vHH = vBB+vGG;
      mHH=vHH.M();
      ptHH=vHH.Pt();
    }
    else {
      mHH=-999;
      ptHH=-999;
    }

    if (vRecoJet_tt1.Pt()>0 && vRecoJet_tt2.Pt()>0) {
      LorentzVector vJJ = vRecoJet_tt1+vRecoJet_tt2;
      mJJ_tt=vJJ.M();
      dEta_tt=vRecoJet_tt1.Eta()-vRecoJet_tt2.Eta();
      
      for (Int_t iJet=0; iJet<branchJet->GetEntries(); iJet++) { // reconstructed jet loop                                                                                                  
	jet = (Jet*) branchJet->At(iJet);
	
	if (fabs(jet->Eta)>4.7) continue;
	if (jet->PT<30) continue;
	
	if ( (vRecoJet_tt1.Eta() > vRecoJet_tt2.Eta()) && (jet->Eta > vRecoJet_tt2.Eta()) && (vRecoJet_tt1.Eta() > jet->Eta) ) {
	  nCentral++;
	}
	else if ( (vRecoJet_tt2.Eta() > vRecoJet_tt1.Eta()) && (jet->Eta > vRecoJet_tt1.Eta()) && (vRecoJet_tt2.Eta() > jet->Eta) ) {
	  nCentral++;
	}
      }	
    }
    else {
      mJJ_tt=-999;
      dEta_tt=0;
    }

    if (vRecoJet_6j1.Pt()>0 && vRecoJet_6j2.Pt()>0) {
      LorentzVector vJJ = vRecoJet_6j1+vRecoJet_6j2;
      mJJ_6j=vJJ.M();
      dEta_6j=vRecoJet_6j1.Eta()-vRecoJet_6j2.Eta();
      
      if (vRecoJet_6j1.Eta() > vRecoJet_6j2.Eta()) {
	if (vRecoB1.Eta() > vRecoJet_6j2.Eta() && vRecoJet_6j1.Eta() > vRecoB1.Eta()) centB++;
	if (vRecoB2.Eta() > vRecoJet_6j2.Eta() && vRecoJet_6j1.Eta() > vRecoB2.Eta()) centB++;
	if (vRecoB3.Eta() > vRecoJet_6j2.Eta() && vRecoJet_6j1.Eta() > vRecoB3.Eta()) centB++;
	if (vRecoB4.Eta() > vRecoJet_6j2.Eta() && vRecoJet_6j1.Eta() > vRecoB4.Eta()) centB++;	
      }	
      else if (vRecoJet_6j2.Eta() > vRecoJet_6j1.Eta()) {
	if (vRecoB1.Eta() > vRecoJet_6j1.Eta() && vRecoJet_6j2.Eta() > vRecoB1.Eta()) centB++;
	if (vRecoB2.Eta() > vRecoJet_6j1.Eta() && vRecoJet_6j2.Eta() > vRecoB2.Eta()) centB++;
	if (vRecoB3.Eta() > vRecoJet_6j1.Eta() && vRecoJet_6j2.Eta() > vRecoB3.Eta()) centB++;
	if (vRecoB4.Eta() > vRecoJet_6j1.Eta() && vRecoJet_6j2.Eta() > vRecoB4.Eta()) centB++;	
      }
    }
    else {
      mJJ_6j=-999;
      dEta_6j=0;
    }

    // ********************
    // MT2 CALC
    // ********************

    missET = (MissingET*) branchMET->At(0);

    met=missET->MET;
    metPhi=missET->Phi;

    if (branchPuppiMET) {

      missET = (MissingET*) branchPuppiMET->At(0);
      
      ppMet=missET->MET;
      ppMetPhi=missET->Phi;
    }
    else {
      ppMet=-999;
      ppMetPhi=-999;
    }

    if ( sRecoTau1->Pt()>0 && sRecoTau2->Pt()>0 && sRecoB1->Pt()>0 && sRecoB2->Pt()>0) {

      tau1.SetMagPhi(sRecoTau1->Pt(), sRecoTau1->Phi());
      tau2.SetMagPhi(sRecoTau2->Pt(), sRecoTau2->Phi());
      mTau1=sRecoTau1->M();
      mTau2=sRecoTau2->M();
      
      b1.SetMagPhi(sRecoB1->Pt(), sRecoB1->Phi());
      b2.SetMagPhi(sRecoB2->Pt(), sRecoB2->Phi());
      mB1=sRecoB1->M();
      mB2=sRecoB2->M();
      
      mpt.SetMagPhi(met, metPhi);
      
      TVector2 sumPt = tau1+tau2+mpt;
      
      smT2 calcmt2 = smT2();
      calcmt2.SetB1(b1);
      calcmt2.SetB2(b2);
      calcmt2.SetMPT(sumPt);
      calcmt2.SetMB1(mB1);
      calcmt2.SetMB2(mB2);
      calcmt2.SetMT1(mTau1);
      calcmt2.SetMT2(mTau2);

      TVector2 c1=sumPt;
      TVector2 c2=sumPt-c1;
      
      ROOT::Math::Functor f(calcmt2,2);
      double step[2] = {0.1, 0.1};
      double variable[2] = { 0.5*c1.Mod(), 0.0 };
      
      min->SetFunction(f);
      min->SetLimitedVariable(0,"cT",variable[0], step[0], 0.0, sumPt.Mod());
      min->SetLimitedVariable(1,"cPhi",variable[1], step[1], 0.0, TMath::Pi());
      
      min->Minimize();
      mt2 = min->MinValue();

      if (branchPuppiMET) {
	
	ppmpt.SetMagPhi(ppMet, ppMetPhi);
	TVector2 ppSumPt = tau1+tau2+ppmpt;
	
	calcmt2.SetMPT(ppSumPt);
	c1=ppSumPt;
	c2=sumPt-c1;
	
	ROOT::Math::Functor f2(calcmt2,2);
	step[0] = 0.1; step[1] = 0.1;
	variable[0] = 0.5*c1.Mod(); variable[1] = 0.0;
	
	min->SetFunction(f2);
	min->SetLimitedVariable(0,"cT",variable[0], step[0], 0.0, ppSumPt.Mod());
	min->SetLimitedVariable(1,"cPhi",variable[1], step[1], 0.0, TMath::Pi());
	
	min->Minimize();
	ppMt2 = min->MinValue();
      }
      else ppMt2 = -999;
    }    
    else { 
      mt2 = -999;
      ppMt2 = -999;
    }
    
    if (ptB1>0 && ptB2>0 && ptTau1>0 && ptTau2>0) { isBBTT=1; }
    if (ptB1>0 && ptB2>0 && ptB3>0 && ptB4>0) { isBBBB=1; }
    if (ptB1>0 && ptB2>0 && ptG1>0 && ptG2>0) { isBBGG=1; }
    if (ptJet_tt1>0 && ptJet_tt2>0 && ptTau1>0 && ptTau2>0) { isVBFTT=1; }
    if (ptB1>0 && ptB2>0 && ptB3>0 && ptB4>0 && ptJet_6j1>0 && ptJet_6j2>0) { isVBF4B=1; }

    if (isBBTT==0 && isBBBB==0 && isBBGG==0 && isVBFTT==0 && isVBF4B==0) continue;
    
    // ********************
    // GEN PARTICLES
    // ********************

    Int_t nH=0, nW=0, nZ=0, nT=0, nB=0, nQ=0, nG=0, nP=0, nL=0;
    for (Int_t iParticle=0; iParticle<branchParticle->GetEntries(); iParticle++) { // generator particle loop
      genParticle = (GenParticle*) branchParticle->At(iParticle);
      
      Int_t pid=fabs(genParticle->PID);
      if (pid==2212) continue;
      if (pid==25) {
	nH++;
	if (iH1==-1) {
	  iH1 = iParticle;
	  genH1 = (GenParticle*) branchParticle->At(iH1);
	}
	else if (iH2==-1) {
	  iH2 = iParticle;
	  genH2 = (GenParticle*) branchParticle->At(iH2);
	}
      }
      else if (pid==23) nZ++;
      else if (pid==24) nW++;
      else if (pid==6) nT++;
      else if (pid==5) nB++;
      else if (pid<5)  nQ++;
      else if (pid==21) nG++;
      else if (pid==22) nP++;
      else if (pid<17 && pid>10) nL++;

      if ( fabs(genParticle->PID) == TAU_ID_CODE ) { // tau switch
	if ( deltaR(genParticle->Eta, vRecoTau1.Eta(), genParticle->Phi, vRecoTau1.Phi()) < MAX_MATCH_DIST ) {
	  iGenTau1 = iParticle;
	  genTau1 = (GenParticle*) branchParticle->At(iGenTau1);
	}
	else if ( deltaR(genParticle->Eta, vRecoTau2.Eta(), genParticle->Phi, vRecoTau2.Phi()) < MAX_MATCH_DIST ) {
	  iGenTau2 = iParticle;
	  genTau2 = (GenParticle*) branchParticle->At(iGenTau2);
	}
	else if ( (vRecoTau1.Pt()<0) && (iGenTau1==-1) ) { 
	  iGenTau1 = iParticle;
	  genTau1 = (GenParticle*) branchParticle->At(iGenTau1);
	}
	else if ( (vRecoTau2.Pt()<0) && (iGenTau2==-1) ) { 
	  iGenTau2 = iParticle;
	  genTau2 = (GenParticle*) branchParticle->At(iGenTau2);
	}
      }

      if ( fabs(genParticle->PID) == G_ID_CODE ) { // gamma switch

	if ( deltaR(genParticle->Eta, vRecoGam1.Eta(), genParticle->Phi, vRecoGam1.Phi()) < MAX_MATCH_DIST ) {
	  iGenGam1 = iParticle;
	  genGam1 = (GenParticle*) branchParticle->At(iGenGam1);
	}
	else if ( deltaR(genParticle->Eta, vRecoGam2.Eta(), genParticle->Phi, vRecoGam2.Phi()) < MAX_MATCH_DIST ) {
	  iGenGam2 = iParticle;
	  genGam2 = (GenParticle*) branchParticle->At(iGenGam2);
	}
	else if ( (vRecoGam1.Pt()<0) && (iGenGam1==-1) ) { 
	  iGenGam1 = iParticle;
	  genGam1 = (GenParticle*) branchParticle->At(iGenGam1);
	}
	else if ( (vRecoGam2.Pt()<0) && (iGenGam2==-1) ) { 
	  iGenGam2 = iParticle;
	  genGam2 = (GenParticle*) branchParticle->At(iGenGam2);
	}
      }
      
      if ( fabs(genParticle->PID) == B_ID_CODE ) { // b switch
	
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
	else if (vRecoB1.Pt()<0) { 
	  iGenB1 = iParticle;
	  genB1 = (GenParticle*) branchParticle->At(iGenB1);
	}
	else if (vRecoB2.Pt()<0) { 
	  iGenB2 = iParticle;
	  genB2 = (GenParticle*) branchParticle->At(iGenB2);
	}
	else if (vRecoB3.Pt()<0) { 
	  iGenB3 = iParticle;
	  genB3 = (GenParticle*) branchParticle->At(iGenB3);
	}
	else if (vRecoB4.Pt()<0) { 
	  iGenB4 = iParticle;
	  genB4 = (GenParticle*) branchParticle->At(iGenB4);
	}
      }
    }

    LorentzVector vGenTau1(0,0,0,0);
    if (genTau1) {
      vGenTau1.SetPt(genTau1->PT);
      vGenTau1.SetEta(genTau1->Eta);
      vGenTau1.SetPhi(genTau1->Phi);
      vGenTau1.SetM(genTau1->Mass);
      sGenTau1 = &vGenTau1;
    }
    else sGenTau1 = &nothing;

    LorentzVector vGenTau2(0,0,0,0);
    if (genTau2) {
      vGenTau2.SetPt(genTau2->PT);
      vGenTau2.SetEta(genTau2->Eta);
      vGenTau2.SetPhi(genTau2->Phi);
      vGenTau2.SetM(genTau2->Mass);
      sGenTau2 = &vGenTau2;
    }
    else sGenTau2 = &nothing;

    LorentzVector vGenGam1(0,0,0,0);
    if (genGam1) {
      vGenGam1.SetPt(genGam1->PT);
      vGenGam1.SetEta(genGam1->Eta);
      vGenGam1.SetPhi(genGam1->Phi);
      vGenGam1.SetM(0);
      sGenGam1 = &vGenGam1;
    }
    else sGenGam1 = &nothing;

    LorentzVector vGenGam2(0,0,0,0);
    if (genGam2) {
      vGenGam2.SetPt(genGam2->PT);
      vGenGam2.SetEta(genGam2->Eta);
      vGenGam2.SetPhi(genGam2->Phi);
      vGenGam2.SetM(0);
      sGenGam2 = &vGenGam2;
    }
    else sGenGam2 = &nothing;

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


    // ********************
    // GEN JETS
    // ********************    

    for (Int_t iJet=0; iJet<branchGenJet->GetEntries(); iJet++) { // generator level jet loop                                                                                      
      genJet = (Jet*) branchGenJet->At(iJet);

      if (deltaR(genJet->Eta, sRecoB1->Eta(), genJet->Phi, sRecoB1->Phi()) < MAX_MATCH_DIST) {
        iGenJetB1=iJet;
        genJetB1 = (Jet*) branchGenJet->At(iGenJetB1);
      }
      else if (deltaR(genJet->Eta, sRecoB2->Eta(), genJet->Phi, sRecoB2->Phi()) < MAX_MATCH_DIST) {
        iGenJetB2=iJet;
        genJetB2 = (Jet*) branchGenJet->At(iGenJetB2);
      }
      else if (deltaR(genJet->Eta, sRecoB3->Eta(), genJet->Phi, sRecoB3->Phi()) < MAX_MATCH_DIST) {
        iGenJetB3=iJet;
        genJetB3 = (Jet*) branchGenJet->At(iGenJetB3);
      }
      else if (deltaR(genJet->Eta, sRecoB4->Eta(), genJet->Phi, sRecoB4->Phi()) < MAX_MATCH_DIST) {
        iGenJetB4=iJet;
        genJetB4 = (Jet*) branchGenJet->At(iGenJetB4);
      }
      else if (deltaR(genJet->Eta, sRecoTau1->Eta(), genJet->Phi, sRecoTau1->Phi()) < MAX_MATCH_DIST) {
        iGenJetTau1=iJet;
        genJetTau1 = (Jet*) branchGenJet->At(iGenJetTau1);
      }
      else if (deltaR(genJet->Eta, sRecoTau2->Eta(), genJet->Phi, sRecoTau2->Phi()) < MAX_MATCH_DIST) {
        iGenJetTau2=iJet;
        genJetTau2 = (Jet*) branchGenJet->At(iGenJetTau2);
      }
      else if (deltaR(genJet->Eta, sRecoJet_tt1->Eta(), genJet->Phi, sRecoJet_tt1->Phi()) < MAX_MATCH_DIST) {
        iGenJet_tt1=iJet;
        genJet_tt1 = (Jet*) branchGenJet->At(iGenJet_tt1);
      }
      else if (deltaR(genJet->Eta, sRecoJet_tt2->Eta(), genJet->Phi, sRecoJet_tt2->Phi()) < MAX_MATCH_DIST) {
        iGenJet_tt2=iJet;
        genJet_tt2 = (Jet*) branchGenJet->At(iGenJet_tt2);
      }
      else if (deltaR(genJet->Eta, sRecoJet_6j1->Eta(), genJet->Phi, sRecoJet_6j1->Phi()) < MAX_MATCH_DIST) {
        iGenJet_6j1=iJet;
        genJet_6j1 = (Jet*) branchGenJet->At(iGenJet_6j1);
      }
      else if (deltaR(genJet->Eta, sRecoJet_6j2->Eta(), genJet->Phi, sRecoJet_6j2->Phi()) < MAX_MATCH_DIST) {
        iGenJet_6j2=iJet;
        genJet_6j2 = (Jet*) branchGenJet->At(iGenJet_6j2);
      }
    }

    LorentzVector vGenJetTau1(0,0,0,0);
    if (genJetTau1) {
      vGenJetTau1.SetPt(genJetTau1->PT);
      vGenJetTau1.SetEta(genJetTau1->Eta);
      vGenJetTau1.SetPhi(genJetTau1->Phi);
      vGenJetTau1.SetM(genJetTau1->Mass);
      sGenJetTau1 = &vGenJetTau1;
    }
    else sGenJetTau1 = &nothing;

    LorentzVector vGenJetTau2(0,0,0,0);
    if (genJetTau2) {
      vGenJetTau2.SetPt(genJetTau2->PT);
      vGenJetTau2.SetEta(genJetTau2->Eta);
      vGenJetTau2.SetPhi(genJetTau2->Phi);
      vGenJetTau2.SetM(genJetTau2->Mass);
      sGenJetTau2 = &vGenJetTau2;
    }
    else sGenJetTau2 = &nothing;
    
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

    LorentzVector vGenJet_tt1(0,0,0,0);
    if (genJet_tt1) {
      vGenJet_tt1.SetPt(genJet_tt1->PT);
      vGenJet_tt1.SetEta(genJet_tt1->Eta);
      vGenJet_tt1.SetPhi(genJet_tt1->Phi);
      vGenJet_tt1.SetM(genJet_tt1->Mass);
      sGenJet_tt1 = &vGenJet_tt1;
    }
    else sGenJet_tt1 = &nothing;

    LorentzVector vGenJet_tt2(0,0,0,0);
    if (genJet_tt2) {
      vGenJet_tt2.SetPt(genJet_tt2->PT);
      vGenJet_tt2.SetEta(genJet_tt2->Eta);
      vGenJet_tt2.SetPhi(genJet_tt2->Phi);
      vGenJet_tt2.SetM(genJet_tt2->Mass);
      sGenJet_tt2 = &vGenJet_tt2;
    }
    else sGenJet_tt2 = &nothing;

    LorentzVector vGenJet_6j1(0,0,0,0);
    if (genJet_6j1) {
      vGenJet_6j1.SetPt(genJet_6j1->PT);
      vGenJet_6j1.SetEta(genJet_6j1->Eta);
      vGenJet_6j1.SetPhi(genJet_6j1->Phi);
      vGenJet_6j1.SetM(genJet_6j1->Mass);
      sGenJet_6j1 = &vGenJet_6j1;
    }
    else sGenJet_6j1 = &nothing;

    LorentzVector vGenJet_6j2(0,0,0,0);
    if (genJet_6j2) {
      vGenJet_6j2.SetPt(genJet_6j2->PT);
      vGenJet_6j2.SetEta(genJet_6j2->Eta);
      vGenJet_6j2.SetPhi(genJet_6j2->Phi);
      vGenJet_6j2.SetM(genJet_6j2->Mass);
      sGenJet_6j2 = &vGenJet_6j2;
    }
    else sGenJet_6j2 = &nothing;

    if ( (isBBBB==1 || isVBF4B==1) && (vGenH1.Pt()>0 && vGenH2.Pt()>0) ) {

      LorentzVector vTestBB1=vGenB1+vGenB2;
      LorentzVector vTestBB2=vGenB3+vGenB4;
      
      if ( (deltaR(vTestBB1.Eta(), vGenH1.Eta(), vTestBB1.Phi(), vGenH1.Phi()) < MAX_MATCH_DIST) && (deltaR(vTestBB2.Eta(), vGenH2.Eta(), vTestBB2.Phi(), vGenH2.Phi()) < MAX_MATCH_DIST) ) {
	cout << "H matched " << endl;
	iHmatch1=1;
	iHmatch2=1;
	iHmatch3=2;
	iHmatch4=2;
      }
      else if ( (deltaR(vTestBB2.Eta(), vGenH1.Eta(), vTestBB2.Phi(), vGenH1.Phi()) < MAX_MATCH_DIST) && (deltaR(vTestBB1.Eta(), vGenH2.Eta(), vTestBB1.Phi(), vGenH2.Phi()) < MAX_MATCH_DIST) ) {
	cout << "H matched " << endl;
	iHmatch1=2;
	iHmatch2=2;
	iHmatch3=1;
	iHmatch4=1;
      }
      
      vTestBB1=vGenB1+vGenB3;
      vTestBB2=vGenB2+vGenB4;
      
      if ( (deltaR(vTestBB1.Eta(), vGenH1.Eta(), vTestBB1.Phi(), vGenH1.Phi()) < MAX_MATCH_DIST) && (deltaR(vTestBB2.Eta(), vGenH2.Eta(), vTestBB2.Phi(), vGenH2.Phi()) < MAX_MATCH_DIST) ) {
	cout << "H matched " << endl;
	iHmatch1=1;
	iHmatch2=2;
	iHmatch3=1;
	iHmatch4=2;
      }
      else if ( (deltaR(vTestBB2.Eta(), vGenH1.Eta(), vTestBB2.Phi(), vGenH1.Phi()) < MAX_MATCH_DIST) && (deltaR(vTestBB1.Eta(), vGenH2.Eta(), vTestBB1.Phi(), vGenH2.Phi()) < MAX_MATCH_DIST) ) {
	cout << "H matched " << endl;
	iHmatch1=2;
	iHmatch2=1;
	iHmatch3=2;
	iHmatch4=1;
      }
      
      vTestBB1=vGenB1+vGenB4;
      vTestBB2=vGenB3+vGenB2;

      if ( (deltaR(vTestBB1.Eta(), vGenH1.Eta(), vTestBB1.Phi(), vGenH1.Phi()) < MAX_MATCH_DIST) && (deltaR(vTestBB2.Eta(), vGenH2.Eta(), vTestBB2.Phi(), vGenH2.Phi()) < MAX_MATCH_DIST) ) {
	cout << "H matched " << endl;
	iHmatch1=1;
	iHmatch2=2;
	iHmatch3=2;
	iHmatch4=1;
      }
      else if ( (deltaR(vTestBB2.Eta(), vGenH1.Eta(), vTestBB2.Phi(), vGenH1.Phi()) < MAX_MATCH_DIST) && (deltaR(vTestBB1.Eta(), vGenH2.Eta(), vTestBB1.Phi(), vGenH2.Phi()) < MAX_MATCH_DIST) ) {
	cout << "H matched " << endl;
	iHmatch1=2;
	iHmatch2=1;
	iHmatch3=1;
	iHmatch4=2;
      }
    }
    
    if (iHmatch1==-1) {
      iHmatch1=0;
      iHmatch2=0;
      iHmatch3=0;
      iHmatch4=0;
    }

    // ********************
    // BKGD SORTING
    // ********************

    if (nH>10) nH=9; if (nW>10) nW=9;
    if (nZ>10) nZ=9; if (nT>10) nT=9;
    if (nB>10) nB=9; if (nQ>10) nQ=9;
    if (nG>10) nG=9; if (nP>10) nP=9;
    if (nL>10) nL=9;

    genInfo=nH+1e1*nW+1e2*nZ+1e3*nT+1e4*nB+1e5*nQ+1e6*nG+1e7*nP+1e8*nL;

    if ( nH>2 ) eventType=HH;
    else if ( nH>0 ) eventType=H;
    else if ( nT==2 ) eventType=TT;
    else if ( nZ>0 && nW==0 && (nG+nQ)>0 ) eventType=ZJET;
    else if ( nW>0 && nZ==0 && (nG+nQ)>0 ) eventType=WJET;
    else if ( nW+nZ+nT+nL>0 ) eventType=EWK;
    else eventType=ETC;

    outTree->Fill();

  } // end event loop

  outFile->Write();
  outFile->Save();

  cout << endl;

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

  if (fabs(eta)<1.5) {
    if ((meanSqDeltaR<MeanSqDeltaRMaxBarrel)&&(betastar<BetaMinBarrel)) {
      return 0;
    }
    else {
      return 1;
    }
  }
  else if (fabs(eta)<4.0) {
    if ((meanSqDeltaR<MeanSqDeltaRMaxEndcap)&&(betastar<BetaMinEndcap)) {
      return 0;
    }
    else {
      return 1;
    }
  }

  return 1;

}
