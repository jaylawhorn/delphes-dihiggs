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
#include "TLorentzVector.h"
#include "modules/Delphes.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "classes/DelphesClasses.h"
#include "mt2.hh"
#include "TauAnalysis/SVFitHelper/interface/TSVfit.h"
#include "TauAnalysis/SVFitHelper/interface/TSVfitter.h"
#include "JEC/JECHelper/interface/jcorr.h"

#endif

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;

Float_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 );

Float_t deltaPhi( const Float_t phi1, const Float_t phi2 );

Int_t puJetID( Float_t eta, Float_t meanSqDeltaR, Float_t beta );

double doJetcorr(mithep::jcorr *corrector,Jet* ijet,double rho_2,double rho_1,double rho_0);

void selection(const TString inputfile="root://eoscms.cern.ch//store/group/upgrade/delphes/dihiggs_signal_bbtt/gFHHTobbtautau_TuneZ2_8TeV-madgraph/files-v2/hhbbtt_77.root",
	       const Float_t xsec=2.92,
	       Int_t sampleNo=100,
	       const TString outputfile="test.root") {
  
  // declare constants
  const Double_t MUON_MASS = 0.105658369;
  const Double_t ELE_MASS  = 0.000511;
  const Double_t TAU_MASS  = 1.77682;

  const Int_t TAU_ID_CODE = 15;
  const Int_t B_ID_CODE = 5;
  const Int_t G_ID_CODE = 21;

  const Float_t MAX_MATCH_DIST = 0.4;

  const Float_t lcov00 = 2500;
  const Float_t lcov10 = 500;
  const Float_t lcov01 = 500;
  const Float_t lcov11 = 2500;
  const Float_t lcov00pp = 225;
  const Float_t lcov10pp = 45;
  const Float_t lcov01pp = 45;
  const Float_t lcov11pp = 225;

  // event categories
  enum { HH=0, H, TT, WJET, ZJET, EWK, ETC };

  // tau decay modes
  enum { hadron=1, electron, muon };

  // setup mt2 minimizer
  ROOT::Math::Minimizer *min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
  min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
  min->SetTolerance(10.0);
  min->SetPrintLevel(0);

  // setup svfit
  mithep::TSVfitter *fitter = new mithep::TSVfitter();

  TVector2 tau1(0,0), tau2(0,0), mpt(0,0), ppmpt(0,0), pumpt(0,0);
  TVector2 b1(0,0), b2(0,0);
  Double_t mt2=0;
  Double_t mt2puppi=0;
  Double_t mt2pileup=0;
  Double_t m_sv=0; 
  Double_t m_svpileup=0; 
  Double_t m_svpuppi=0;

  //setup jet corrections
  
  mithep::jcorr *corrector = new mithep::jcorr;
  char *PATH = getenv("CMSSW_BASE"); assert(PATH);
  TString path(TString::Format("%s/src/JEC/JECHelper/data/JetCorrections_phase1/", PATH));
  corrector->AddJetCorr(path+"Delphes_V2_MC_L1FastJet_AK4PF.txt");
  corrector->AddJetCorr(path+"Delphes_V2_MC_L2Relative_AK4PF.txt");
  corrector->AddJetCorr(path+"Delphes_V2_MC_L3Absolute_AK4PF.txt");
  corrector->setup();

  // read input input file
  TChain chain("Delphes");
  chain.Add(inputfile);
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  if (!(branchJet)) {
    cout << "File broken" << endl;
    return;
  }

  TClonesArray *branchRawJet = treeReader->UseBranch("RawJet");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
  TClonesArray *branchMET =treeReader->UseBranch("MissingET");
  TClonesArray *branchPuppiMET =treeReader->UseBranch("PuppiMissingET");
  TClonesArray *branchPileupMET =treeReader->UseBranch("PileUpJetIDMissingET");
  TClonesArray *branchRho = treeReader->UseBranch("Rho");
  TClonesArray *branchGenJet = treeReader->UseBranch("GenJet");
  TClonesArray *branchIsoTrack = treeReader->UseBranch("IsoTrack");
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchEvent = treeReader->UseBranch("Event");

  //set up loop variables
  GenParticle *genParticle=0;
  Jet *genJet=0;
  Jet *jet=0;
  Rho *rho0=0;
  Rho *rho1=0;
  Rho *rho2=0;
  Electron *ele=0;
  IsoTrack *iso=0;
  Muon *mu=0;
  Photon *gam=0;
  MissingET *missET=0;
  LHEFEvent *event=0;

  // set up storage variables
  // TAUS
  Jet *jetTau1=0, *jetTau2=0; Electron *eleTau=0; Muon *muTau=0;
  IsoTrack *isoTau1=0, *isoTau2=0;
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
  // GEN TAU-JETS
  Jet *genJetTau1=0, *genJetTau2=0;
  // GEN VBF-JETS
  Jet *genJet_tt1=0, *genJet_tt2=0, *genJet_6j1=0, *genJet_6j2=0;

  Int_t iT1=-1, iT2=-1;
  Int_t iI1=-1, iI2=-1;
  Int_t iB1=-1, iB2=-1, iB3=-1, iB4=-1;
  Int_t iG1=-1, iG2=-1;
  Int_t iJ1=-1, iJ2=-1, iJ3=-1, iJ4=-1;
  Int_t iH1=-1, iH2=-1;

  Int_t iGenTau1=-1, iGenTau2=-1;
  Int_t iGenB1=-1,   iGenB2=-1, iGenB3=-1, iGenB4=-1;
  Int_t iGenGam1=-1, iGenGam2=-1;

  Int_t iGenJetTau1=-1, iGenJetTau2=-1;
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
  Float_t pileupMet, pileupMetPhi;

  Int_t nCentral=0, nBtag=0, nJets=0;
  Int_t centB=0;

  Int_t nLep=0;

  Int_t tauCat1=0, tauCat2=0;
  Int_t bTag1=0, bTag2=0, bTag3=0, bTag4=0;
  Int_t jbTag_tt1=0, jbTag_tt2=0, jbTag_6j1=0, jbTag_6j2=0;
  
  Float_t ptTau1, ptTau2, ptTrk1, ptTrk2, ptB1, ptB2, ptB3, ptB4, ptG1, ptG2, ptJet_tt1, ptJet_tt2, ptJet_6j1, ptJet_6j2, tauIso1, tauIso2;
  Float_t etaTau1, etaTau2, etaTrk1, etaTrk2, etaB1, etaB2, etaB3, etaB4, etaG1, etaG2, etaJet_tt1, etaJet_tt2, etaJet_6j1, etaJet_6j2;
  Float_t phiTau1, phiTau2, phiTrk1, phiTrk2, phiB1, phiB2, phiB3, phiB4, phiG1, phiG2, phiJet_tt1, phiJet_tt2, phiJet_6j1, phiJet_6j2;
  Float_t mTau1, mTau2, mTrk1, mTrk2, mB1, mB2, mB3, mB4, eG1, eG2, mJet_tt1, mJet_tt2, mJet_6j1, mJet_6j2;

  Float_t ptTau1_gen, ptTau2_gen, ptB1_gen, ptB2_gen, ptB3_gen, ptB4_gen, ptG1_gen, ptG2_gen, ptJet_tt1_gen, ptJet_tt2_gen, ptJet_6j1_gen, ptJet_6j2_gen;
  Float_t etaTau1_gen, etaTau2_gen, etaB1_gen, etaB2_gen, etaB3_gen, etaB4_gen, etaG1_gen, etaG2_gen, etaJet_tt1_gen, etaJet_tt2_gen, etaJet_6j1_gen, etaJet_6j2_gen;
  Float_t phiTau1_gen, phiTau2_gen, phiB1_gen, phiB2_gen, phiB3_gen, phiB4_gen, phiG1_gen, phiG2_gen, phiJet_tt1_gen, phiJet_tt2_gen, phiJet_6j1_gen, phiJet_6j2_gen;
  Float_t mTau1_gen, mTau2_gen, mB1_gen, mB2_gen, mB3_gen, mB4_gen, eG1_gen, eG2_gen, mJet_tt1_gen, mJet_tt2_gen, mJet_6j1_gen, mJet_6j2_gen;

  Float_t ptH1_gen, etaH1_gen, phiH1_gen, mH1_gen;
  Float_t ptH2_gen, etaH2_gen, phiH2_gen, mH2_gen;

  Float_t ptTau1_genJet, ptTau2_genJet, etaTau1_genJet, etaTau2_genJet, phiTau1_genJet, phiTau2_genJet, mTau1_genJet, mTau2_genJet;

  Float_t ptTT, ptBB1, ptBB2, ptB1B2, ptB1B3, ptB1B4, ptB2B3, ptB2B4, ptB3B4, ptGG, ptJJ_tt, ptJJ_6j, ptHH;
  Float_t etaTT, etaBB1, etaBB2, etaB1B2, etaB1B3, etaB1B4, etaB2B3, etaB2B4, etaB3B4, etaGG, etaJJ_tt, etaJJ_6j, etaHH;
  Float_t phiTT, phiBB1, phiBB2, phiB1B2, phiB1B3, phiB1B4, phiB2B3, phiB2B4, phiB3B4,phiGG, phiJJ_tt, phiJJ_6j, phiHH;
  Float_t mTT, mBB1, mBB2, mB1B2, mB1B3, mB1B4, mB2B3, mB2B4, mB3B4,mGG, mJJ_tt, mJJ_6j, mHH;

  Float_t dEta_tt=0, dEta_6j=0;

  Float_t rho_0=0, rho_1=0, rho_2=0;

  Float_t dEtaBB1=0, dEtaBB2=0, dEtaB1B2=0, dEtaB1B3=0, dEtaB1B4=0, dEtaB2B3=0, dEtaB2B4=0, dEtaB3B4=0, dEtaTT=0, dEtaHH=0;
  Float_t dPhiBB1=0, dPhiBB2=0, dPhiB1B2=0, dPhiB1B3=0, dPhiB1B4=0, dPhiB2B3=0, dPhiB2B4=0, dPhiB3B4=0, dPhiTT=0, dPhiHH=0;
  Float_t dRBB1=0, dRBB2=0, dRB1B2=0, dRB1B3=0, dRB1B4=0, dRB2B3=0, dRB2B4=0, dRB3B4=0, dRTT=0, dRHH=0;

  Float_t mindR4B=0;
  Int_t nBJetsComb=0;


  TFile *outFile = new TFile(outputfile, "RECREATE");

  TTree *infoTree = new TTree("Count", "Count");
  Long64_t n = numberOfEntries;
  infoTree->Branch("n", &n, "n/i");
  infoTree->Fill();

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
  outTree->Branch("pileupmet",      &pileupMet,      "pileupMet/f");    // pileup missing transverse energy
  outTree->Branch("pileupmetPhi",   &pileupMetPhi,   "pileupMetPhi/f"); // pileup missing transverse energy phi
  outTree->Branch("puppiMet",       &ppMet,          "ppMet/f");        // PUPPI missing transverse energy
  outTree->Branch("puppiMetPhi",    &ppMetPhi,       "ppMetPhi/f");     // PUPPI missing transverse energy phi

  outTree->Branch("ptTau1",         &ptTau1,         "ptTau1/f");       // pt(Tau1)
  outTree->Branch("etaTau1",        &etaTau1,        "etaTau1/f");      // eta(Tau1)
  outTree->Branch("phiTau1",        &phiTau1,        "phiTau1/f");      // phi(Tau1)
  outTree->Branch("mTau1",          &mTau1,          "mTau1/f");        // m(Tau1)
  outTree->Branch("tauCat1",        &tauCat1,        "tauCat1/i");      // leading tau final state - jet, muon, electron
  outTree->Branch("tauIso1",        &tauIso1,        "tauIso1/f");      // leading tau final state - jet, muon, electron

  outTree->Branch("ptTau2",         &ptTau2,         "ptTau2/f");       // pt(Tau2)
  outTree->Branch("etaTau2",        &etaTau2,        "etaTau2/f");      // eta(Tau2)
  outTree->Branch("phiTau2",        &phiTau2,        "phiTau2/f");      // phi(Tau2)
  outTree->Branch("mTau2",          &mTau2,          "mTau2/f");        // m(Tau2)
  outTree->Branch("tauCat2",        &tauCat2,        "tauCat2/i");      // second tau final state - jet, muon, electron
  outTree->Branch("tauIso2",        &tauIso2,        "tauIso2/f");      // second tau final state - jet, muon, electron

  outTree->Branch("ptTrk1",         &ptTrk1,         "ptTrk1/f");       // pt(Trk1)
  outTree->Branch("etaTrk1",        &etaTrk1,        "etaTrk1/f");      // eta(Trk1)
  outTree->Branch("phiTrk1",        &phiTrk1,        "phiTrk1/f");      // phi(Trk1)
  outTree->Branch("mTrk1",          &mTrk1,          "mTrk1/f");        // m(Trk1)

  outTree->Branch("ptTrk2",         &ptTrk2,         "ptTrk2/f");       // pt(Trk2)
  outTree->Branch("etaTrk2",        &etaTrk2,        "etaTrk2/f");      // eta(Trk2)
  outTree->Branch("phiTrk2",        &phiTrk2,        "phiTrk2/f");      // phi(Trk2)
  outTree->Branch("mTrk2",          &mTrk2,          "mTrk2/f");        // m(Trk2)

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

  outTree->Branch("ptB2",           &ptB2,           "ptB2/f");         // pt(B2)
  outTree->Branch("etaB2",          &etaB2,          "etaB2/f");        // eta(B2)
  outTree->Branch("phiB2",          &phiB2,          "phiB2/f");        // phi(B2)
  outTree->Branch("mB2",            &mB2,            "mB2/f");          // m(B2)
  outTree->Branch("bTag2",          &bTag2,          "bTag2/i");        // second b-jet tag from delphes

  outTree->Branch("ptB3",           &ptB3,           "ptB3/f");         // pt(B3)
  outTree->Branch("etaB3",          &etaB3,          "etaB3/f");        // eta(B3)
  outTree->Branch("phiB3",          &phiB3,          "phiB3/f");        // phi(B3)
  outTree->Branch("mB3",            &mB3,            "mB3/f");          // m(B3)
  outTree->Branch("bTag3",          &bTag3,          "bTag3/i");        // third b-jet tag from delphes

  outTree->Branch("ptB4",           &ptB4,           "ptB4/f");         // pt(B4)
  outTree->Branch("etaB4",          &etaB4,          "etaB4/f");        // eta(B4)
  outTree->Branch("phiB4",          &phiB4,          "phiB4/f");        // phi(B4)
  outTree->Branch("mB4",            &mB4,            "mB4/f");          // m(B4)
  outTree->Branch("bTag4",          &bTag4,          "bTag4/i");        // fourth b-jet tag from delphes

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

  outTree->Branch("ptTau1_gen",     &ptTau1_gen,     "ptTau1_gen/f");       // gen pt(Tau1)
  outTree->Branch("etaTau1_gen",    &etaTau1_gen,    "etaTau1_gen/f");      // gen eta(Tau1)
  outTree->Branch("phiTau1_gen",    &phiTau1_gen,    "phiTau1_gen/f");      // gen phi(Tau1)
  outTree->Branch("mTau1_gen",      &mTau1_gen,      "mTau1_gen/f");        // gen m(Tau1)

  outTree->Branch("ptTau2_gen",     &ptTau2_gen,     "ptTau2_gen/f");       // gen pt(Tau2)
  outTree->Branch("etaTau2_gen",    &etaTau2_gen,    "etaTau2_gen/f");      // gen eta(Tau2)
  outTree->Branch("phiTau2_gen",    &phiTau2_gen,    "phiTau2_gen/f");      // gen phi(Tau2)
  outTree->Branch("mTau2_gen",      &mTau2_gen,      "mTau2_gen/f");        // gen m(Tau2)

  outTree->Branch("ptTau1_genJet",  &ptTau1_genJet,  "ptTau1_genJet/f");    // gen pt(Tau1)
  outTree->Branch("etaTau1_genJet", &etaTau1_genJet, "etaTau1_genJet/f");   // gen eta(Tau1)
  outTree->Branch("phiTau1_genJet", &phiTau1_genJet, "phiTau1_genJet/f");   // gen phi(Tau1)
  outTree->Branch("mTau1_genJet",   &mTau1_genJet,   "mTau1_genJet/f");     // gen m(Tau1)

  outTree->Branch("ptTau2_genJet",  &ptTau2_genJet,  "ptTau2_genJet/f");    // gen pt(Tau2)
  outTree->Branch("etaTau2_genJet", &etaTau2_genJet, "etaTau2_genJet/f");   // gen eta(Tau2)
  outTree->Branch("phiTau2_genJet", &phiTau2_genJet, "phiTau2_genJet/f");   // gen phi(Tau2)
  outTree->Branch("mTau2_genJet",   &mTau2_genJet,   "mTau2_genJet/f");     // gen m(Tau2)

  outTree->Branch("ptG1_gen",       &ptG1_gen,       "ptG1_gen/f");         // gen pt(Gam1)
  outTree->Branch("etaG1_gen",      &etaG1_gen,      "etaG1_gen/f");        // gen eta(Gam1)
  outTree->Branch("phiG1_gen",      &phiG1_gen,      "phiG1_gen/f");        // gen phi(Gam1)
  outTree->Branch("eG1_gen",        &eG1_gen,        "eG1_gen/f");          // gen e(Gam1)

  outTree->Branch("ptG2_gen",       &ptG2_gen,       "ptG2_gen/f");         // gen pt(Gam2)
  outTree->Branch("etaG2_gen",      &etaG2_gen,      "etaG2_gen/f");        // gen eta(Gam2)
  outTree->Branch("phiG2_gen",      &phiG2_gen,      "phiG2_gen/f");        // gen phi(Gam2)
  outTree->Branch("eG2_gen",        &eG2_gen,        "eG2_gen/f");          // gen e(Gam2)

  outTree->Branch("ptB1_gen",       &ptB1_gen,       "ptB1_gen/f");         // gen pt(B1)
  outTree->Branch("etaB1_gen",      &etaB1_gen,      "etaB1_gen/f");        // gen eta(B1)
  outTree->Branch("phiB1_gen",      &phiB1_gen,      "phiB1_gen/f");        // gen phi(B1)
  outTree->Branch("mB1_gen",        &mB1_gen,        "mB1_gen/f");          // gen m(B1)
  outTree->Branch("iHmatch1",       &iHmatch1,       "iHmatch1/i");         // if 4 b's, which of two higgs matched

  outTree->Branch("ptB2_gen",       &ptB2_gen,       "ptB2_gen/f");         // gen pt(B2)
  outTree->Branch("etaB2_gen",      &etaB2_gen,      "etaB2_gen/f");        // gen eta(B2)
  outTree->Branch("phiB2_gen",      &phiB2_gen,      "phiB2_gen/f");        // gen phi(B2)
  outTree->Branch("mB2_gen",        &mB2_gen,        "mB2_gen/f");          // gen m(B2)
  outTree->Branch("iHmatch2",       &iHmatch2,       "iHmatch2/i");         // if 4 b's, which of two higgs matched

  outTree->Branch("ptB3_gen",       &ptB3_gen,       "ptB3_gen/f");         // gen pt(B3)
  outTree->Branch("etaB3_gen",      &etaB3_gen,      "etaB3_gen/f");        // gen eta(B3)
  outTree->Branch("phiB3_gen",      &phiB3_gen,      "phiB3_gen/f");        // gen phi(B3)
  outTree->Branch("mB3_gen",        &mB3_gen,        "mB3_gen/f");          // gen m(B3)
  outTree->Branch("iHmatch3",       &iHmatch3,       "iHmatch3/i");         // if 4 b's, which of two higgs matched

  outTree->Branch("ptB4_gen",       &ptB4_gen,       "ptB4_gen/f");         // gen pt(B4)
  outTree->Branch("etaB4_gen",      &etaB4_gen,      "etaB4_gen/f");        // gen eta(B4)
  outTree->Branch("phiB4_gen",      &phiB4_gen,      "phiB4_gen/f");        // gen phi(B4)
  outTree->Branch("mB4_gen",        &mB4_gen,        "mB4_gen/f");          // gen m(B4)
  outTree->Branch("iHmatch4",       &iHmatch4,       "iHmatch4/i");         // if 4 b's, which of two higgs matched

  outTree->Branch("ptH1_gen",       &ptH1_gen,       "ptH1_gen/f");         // gen pt(H1)
  outTree->Branch("etaH1_gen",      &etaH1_gen,      "etaH1_gen/f");        // gen eta(H1)
  outTree->Branch("phiH1_gen",      &phiH1_gen,      "phiH1_gen/f");        // gen phi(H1)
  outTree->Branch("mH1_gen",        &mH1_gen,        "mH1_gen/f");          // gen m(H1)

  outTree->Branch("ptH2_gen",       &ptH2_gen,       "ptH2_gen/f");         // gen pt(H2)
  outTree->Branch("etaH2_gen",      &etaH2_gen,      "etaH2_gen/f");        // gen eta(H2)
  outTree->Branch("phiH2_gen",      &phiH2_gen,      "phiH2_gen/f");        // gen phi(H2)
  outTree->Branch("mH2_gen",        &mH2_gen,        "mH2_gen/f");          // gen m(H2)

  outTree->Branch("ptJet_tt1_gen",  &ptJet_tt1_gen,  "ptJet_tt1_gen/f");    // gen jet pt(Jet1)
  outTree->Branch("etaJet_tt1_gen", &etaJet_tt1_gen, "etaJet_tt1_gen/f");   // gen jet eta(Jet1)
  outTree->Branch("phiJet_tt1_gen", &phiJet_tt1_gen, "phiJet_tt1_gen/f");   // gen jet phi(Jet1)
  outTree->Branch("mJet_tt1_gen",   &mJet_tt1_gen,   "mJet_tt1_gen/f");     // gen jet m(Jet1)

  outTree->Branch("ptJet_tt2_gen",  &ptJet_tt2_gen,  "ptJet_tt2_gen/f");    // gen jet pt(Jet_tt2)
  outTree->Branch("etaJet_tt2_gen", &etaJet_tt2_gen, "etaJet_tt2_gen/f");   // gen jet eta(Jet_tt2)
  outTree->Branch("phiJet_tt2_gen", &phiJet_tt2_gen, "phiJet_tt2_gen/f");   // gen jet phi(Jet_tt2)
  outTree->Branch("mJet_tt2_gen",   &mJet_tt2_gen,   "mJet_tt2_gen/f");     // gen jet m(Jet_tt2)

  outTree->Branch("ptJet_6j1_gen",  &ptJet_6j1_gen,  "ptJet_6j1_gen/f");    // gen jet pt(Jet1)
  outTree->Branch("etaJet_6j1_gen", &etaJet_6j1_gen, "etaJet_6j1_gen/f");   // gen jet eta(Jet1)
  outTree->Branch("phiJet_6j1_gen", &phiJet_6j1_gen, "phiJet_6j1_gen/f");   // gen jet phi(Jet1)
  outTree->Branch("mJet_6j1_gen",   &mJet_6j1_gen,   "mJet_6j1_gen/f");     // gen jet m(Jet1)

  outTree->Branch("ptJet_6j2_gen",  &ptJet_6j2_gen,  "ptJet_6j2_gen/f");    // gen jet pt(Jet2)
  outTree->Branch("etaJet_6j2_gen", &etaJet_6j2_gen, "etaJet_6j2_gen/f");   // gen jet eta(Jet2)
  outTree->Branch("phiJet_6j2_gen", &phiJet_6j2_gen, "phiJet_6j2_gen/f");   // gen jet phi(Jet2)
  outTree->Branch("mJet_6j2_gen",   &mJet_6j2_gen,   "mJet_6j2_gen/f");     // gen jet m(Jet2)

  outTree->Branch("ptTT",           &ptTT,           "ptTT/f");         // pt(TT)
  outTree->Branch("etaTT",          &etaTT,          "etaTT/f");        // eta(TT)
  outTree->Branch("phiTT",          &phiTT,          "phiTT/f");        // phi(TT)
  outTree->Branch("mTT",            &mTT,            "mTT/f");          // m(TT)

  outTree->Branch("ptBB1",          &ptBB1,          "ptBB1/f");        // pt(BB1)
  outTree->Branch("etaBB1",         &etaBB1,         "etaBB1/f");       // eta(BB1)
  outTree->Branch("phiBB1",         &phiBB1,         "phiBB1/f");       // phi(BB1)
  outTree->Branch("mBB1",           &mBB1,           "mBB1/f");         // m(BB1)

  outTree->Branch("ptBB2",          &ptBB2,          "ptBB2/f");        // pt(BB2)
  outTree->Branch("etaBB2",         &etaBB2,         "etaBB2/f");       // eta(BB2)
  outTree->Branch("phiBB2",         &phiBB2,         "phiBB2/f");       // phi(BB2)
  outTree->Branch("mBB2",           &mBB2,           "mBB2/f");         // m(BB2)

  outTree->Branch("ptB1B2",          &ptB1B2,          "ptB1B2/f");        // pt(B1B2)
  outTree->Branch("etaB1B2",         &etaB1B2,         "etaB1B2/f");       // eta(B1B2)
  outTree->Branch("phiB1B2",         &phiB1B2,         "phiB1B2/f");       // phi(B1B2)
  outTree->Branch("mB1B2",           &mB1B2,           "mB1B2/f");         // m(B1B2)

  outTree->Branch("ptB1B3",          &ptB1B3,          "ptB1B3/f");        // pt(B1B3)
  outTree->Branch("etaB1B3",         &etaB1B3,         "etaB1B3/f");       // eta(B1B3)
  outTree->Branch("phiB1B3",         &phiB1B3,         "phiB1B3/f");       // phi(B1B3)
  outTree->Branch("mB1B3",           &mB1B3,           "mB1B3/f");         // m(B1B3)

  outTree->Branch("ptB1B4",          &ptB1B4,          "ptB1B4/f");        // pt(B1B4)
  outTree->Branch("etaB1B4",         &etaB1B4,         "etaB1B4/f");       // eta(B1B4)
  outTree->Branch("phiB1B4",         &phiB1B4,         "phiB1B4/f");       // phi(B1B4)
  outTree->Branch("mB1B4",           &mB1B4,           "mB1B4/f");         // m(B1B4)

  outTree->Branch("ptB2B3",          &ptB2B3,          "ptB2B3/f");        // pt(B2B3)
  outTree->Branch("etaB2B3",         &etaB2B3,         "etaB2B3/f");       // eta(B2B3)
  outTree->Branch("phiB2B3",         &phiB2B3,         "phiB2B3/f");       // phi(B2B3)
  outTree->Branch("mB2B3",           &mB2B3,           "mB2B3/f");         // m(B2B3)

  outTree->Branch("ptB2B4",          &ptB2B4,          "ptB2B4/f");        // pt(B2B4)
  outTree->Branch("etaB2B4",         &etaB2B4,         "etaB2B4/f");       // eta(B2B4)
  outTree->Branch("phiB2B4",         &phiB2B4,         "phiB2B4/f");       // phi(B2B4)
  outTree->Branch("mB2B4",           &mB2B4,           "mB2B4/f");         // m(B2B4)

  outTree->Branch("ptB3B4",          &ptB3B4,          "ptB3B4/f");        // pt(B3B4)
  outTree->Branch("etaB3B4",         &etaB3B4,         "etaB3B4/f");       // eta(B3B4)
  outTree->Branch("phiB3B4",         &phiB3B4,         "phiB3B4/f");       // phi(B3B4)
  outTree->Branch("mB3B4",           &mB3B4,           "mB3B4/f");         // m(B3B4)

  outTree->Branch("ptGG",           &ptGG,           "ptGG/f");         // pt(GG)
  outTree->Branch("etaGG",          &etaGG,          "etaGG/f");        // eta(GG)
  outTree->Branch("phiGG",          &phiGG,          "phiGG/f");        // phi(GG)
  outTree->Branch("mGG",            &mGG,            "mGG/f");          // m(GG)

  outTree->Branch("ptJJ_tt",        &ptJJ_tt,        "ptJJ_tt/f");      // pt(JJ_tt)
  outTree->Branch("etaJJ_tt",       &etaJJ_tt,       "etaJJ_tt/f");     // eta(JJ_tt)
  outTree->Branch("phiJJ_tt",       &phiJJ_tt,       "phiJJ_tt/f");     // phi(JJ_tt)
  outTree->Branch("mJJ_tt",         &mJJ_tt,         "mJJ_tt/f");       // m(JJ_tt)

  outTree->Branch("ptJJ_6j",        &ptJJ_6j,        "ptJJ_6j/f");      // pt(JJ_6j)
  outTree->Branch("etaJJ_6j",       &etaJJ_6j,       "etaJJ_6j/f");     // eta(JJ_6j)
  outTree->Branch("phiJJ_6j",       &phiJJ_6j,       "phiJJ_6j/f");     // phi(JJ_6j)
  outTree->Branch("mJJ_6j",         &mJJ_6j,         "mJJ_6j/f");       // m(JJ_6j)

  outTree->Branch("ptHH",           &ptHH,           "ptHH/f");         // pt(HH)
  outTree->Branch("etaHH",          &etaHH,          "etaHH/f");        // eta(HH)
  outTree->Branch("phiHH",          &phiHH,          "phiHH/f");        // phi(HH)
  outTree->Branch("mHH",            &mHH,            "mHH/f");          // m(HH)

  outTree->Branch("mt2",            &mt2,            "mt2/D");          // "stransverse mass" (HH)
  outTree->Branch("mt2puppi",       &mt2puppi,       "mt2puppi/D");     // PUPPI "stransverse mass" (HH)
  outTree->Branch("mt2pileup",      &mt2pileup,      "mt2pileup/D");    // PUJetID "stransverse mass" (HH)
  
  outTree->Branch("m_sv",           &m_sv,           "m_sv/D");         // "SVFit mass estimate" 
  outTree->Branch("m_svpileup",     &m_svpileup,     "m_svpileup/D");   // "SVFit mass estimate with pileup jet ID MET"
  outTree->Branch("m_svpuppi",      &m_svpuppi,      "m_svpuppi/D");    // "SVFit mass estimate with pileup jet ID MET"

  outTree->Branch("nBtag",          &nBtag,          "nBtag/i");        // number of b-tagged jets (VBF)   
  outTree->Branch("nCentral",       &nCentral,       "nCentral/i");     // number of central jets (VBF)
  outTree->Branch("centB",          &centB,          "centB/i");        // how many b's between VBF jets?
  outTree->Branch("nJets",          &nJets,          "nJets/i");        // number of jets
  outTree->Branch("nLep",           &nLep,           "nLep/i");        // number of leptons
  outTree->Branch("dEta_tt",        &dEta_tt,        "dEta_tt/f");      // delta Eta (VBF)
  outTree->Branch("dEta_6j",        &dEta_6j,        "dEta_6j/f");      // delta Eta (VBF)
  outTree->Branch("rho_0",          &rho_0,          "rho_0/f");        // central rho, 0-2.5
  outTree->Branch("rho_1",          &rho_1,          "rho_1/f");        // central rho, 2.5-4
  outTree->Branch("rho_2",          &rho_2,          "rho_2/f");        // central rho, 4-5

  outTree->Branch("dEtaBB1",         &dEtaBB1,        "dEtaBB1/f");
  outTree->Branch("dPhiBB1",         &dPhiBB1,        "dPhiBB1/f");
  outTree->Branch("dRBB1",           &dRBB1,          "dRBB1/f");

  outTree->Branch("dEtaBB2",         &dEtaBB2,        "dEtaBB2/f");
  outTree->Branch("dPhiBB2",         &dPhiBB2,        "dPhiBB2/f");
  outTree->Branch("dRBB2",           &dRBB2,          "dRBB2/f");

  outTree->Branch("dEtaB1B2",         &dEtaB1B2,        "dEtaB1B2/f");
  outTree->Branch("dPhiB1B2",         &dPhiB1B2,        "dPhiB1B2/f");
  outTree->Branch("dRB1B2",           &dRB1B2,          "dRB1B2/f");

  outTree->Branch("dEtaB1B3",         &dEtaB1B3,        "dEtaB1B3/f");
  outTree->Branch("dPhiB1B3",         &dPhiB1B3,        "dPhiB1B3/f");
  outTree->Branch("dRB1B3",           &dRB1B3,          "dRB1B3/f");

  outTree->Branch("dEtaB1B4",         &dEtaB1B4,        "dEtaB1B4/f");
  outTree->Branch("dPhiB1B4",         &dPhiB1B4,        "dPhiB1B4/f");
  outTree->Branch("dRB1B4",           &dRB1B4,          "dRB1B4/f");

  outTree->Branch("dEtaB2B3",         &dEtaB2B3,        "dEtaB2B3/f");
  outTree->Branch("dPhiB2B3",         &dPhiB2B3,        "dPhiB2B3/f");
  outTree->Branch("dRB2B3",           &dRB2B3,          "dRB2B3/f");

  outTree->Branch("dEtaB2B4",         &dEtaB2B4,        "dEtaB2B4/f");
  outTree->Branch("dPhiB2B4",         &dPhiB2B4,        "dPhiB2B4/f");
  outTree->Branch("dRB2B4",           &dRB2B4,          "dRB2B4/f");

  outTree->Branch("dEtaB3B4",         &dEtaB3B4,        "dEtaB3B4/f");
  outTree->Branch("dPhiB3B4",         &dPhiB3B4,        "dPhiB3B4/f");
  outTree->Branch("dRB3B4",           &dRB3B4,          "dRB3B4/f");

  outTree->Branch("mindR4B",           &mindR4B,          "mindR4B/f");
  outTree->Branch("nBJetsComb",        &nBJetsComb,       "nBJetsComb/f");

  outTree->Branch("dEtaTT",          &dEtaTT,         "dEtaTT/f");
  outTree->Branch("dPhiTT",          &dPhiTT,         "dPhiTT/f");
  outTree->Branch("dRTT",            &dRTT,           "dRTT/f");

  outTree->Branch("dEtaHH",          &dEtaHH,         "dEtaHH/f");
  outTree->Branch("dPhiHH",          &dPhiHH,         "dPhiHH/f");
  outTree->Branch("dRHH",            &dRHH,           "dRHH/f");

  //for (Int_t iEntry=0; iEntry<10; iEntry++) { // entry loop
  for (Int_t iEntry=0; iEntry<numberOfEntries; iEntry++) { // entry loop
    treeReader->ReadEntry(iEntry);
    //cout << ".";
    //cout << iEntry << endl;

    // ********************
    // RESET
    // ********************

    iB1=-1; iB2=-1; iB3=-1; iB4=-1; iG1=-1; iG2=-1; 
    iT1=-1; iT2=-1; iJ1=-1; iJ2=-1; iJ3=-1; iJ4=-1;
    iH1=-1; iH2=-1; iI1=-1; iI2=-1;

    iGenTau1=-1;    iGenTau2=-1;
    iGenB1=-1;      iGenB2=-1;
    iGenB3=-1;      iGenB4=-1;
    iGenGam1=-1;    iGenGam2=-1;

    iGenJetTau1=-1; iGenJetTau2=-1;
    iGenJet_tt1=-1; iGenJet_tt2=-1;
    iGenJet_6j1=-1; iGenJet_6j2=-1;

    tauCat1=0; tauCat2=0; 
    bTag1=0; bTag2=0; bTag3=0; bTag4=0;
    jbTag_tt1=0; jbTag_tt2=0;
    jbTag_6j1=0; jbTag_6j2=0;
    
    mTT=-999; mBB1=-999; mBB2=-999; mB1B2=-999; mB1B3=-999; mB1B4=-999; mB2B3=-999; mB2B4=-999; mB3B4=-999; mGG=-999; mHH=-999; mJJ_tt=-999; mJJ_6j=-999; 
    ptTT=-999; ptBB1=-999; ptBB2=-999; ptB1B2=-999; ptB1B3=-999; ptB1B4=-999; ptB2B3=-999; ptB2B4=-999; ptB3B4=-999; ptGG=-999; ptHH=-999; ptJJ_tt=-999; ptJJ_6j=-999; 
    etaTT=-999; etaBB1=-999; etaBB2=-999; etaB1B2=-999; etaB1B3=-999; etaB1B4=-999; etaB2B3=-999; etaB2B4=-999; etaB3B4=-999; etaGG=-999; etaHH=-999; etaJJ_tt=-999; etaJJ_6j=-999; 
    phiTT=-999; phiBB1=-999; phiBB2=-999; phiB1B2=-999; phiB1B3=-999; phiB1B4=-999; phiB2B3=-999; phiB2B4=-999; phiB3B4=-999; phiGG=-999; phiHH=-999; phiJJ_tt=-999; phiJJ_6j=-999; 

    ptTau1=-999; etaTau1=-999; phiTau1=-999; mTau1=-999; tauIso1=-999;
    ptTau2=-999; etaTau2=-999; phiTau2=-999; mTau2=-999; tauIso2=-999;
    ptTrk1=-999; etaTrk1=-999; phiTrk1=-999; mTrk1=-999;
    ptTrk2=-999; etaTrk2=-999; phiTrk2=-999; mTrk2=-999;
    ptG1=-999; etaG1=-999; phiG1=-999; eG1=-999;
    ptG2=-999; etaG2=-999; phiG2=-999; eG2=-999;
    ptB1=-999; etaB1=-999; phiB1=-999; mB1=-999;
    ptB2=-999; etaB2=-999; phiB2=-999; mB2=-999;
    ptB3=-999; etaB3=-999; phiB3=-999; mB3=-999;
    ptB4=-999; etaB4=-999; phiB4=-999; mB4=-999;

    ptTau1_gen=-999; etaTau1_gen=-999; phiTau1_gen=-999; mTau1_gen=-999;
    ptTau2_gen=-999; etaTau2_gen=-999; phiTau2_gen=-999; mTau2_gen=-999;
    ptG1_gen=-999; etaG1_gen=-999; phiG1_gen=-999; eG1_gen=-999;
    ptG2_gen=-999; etaG2_gen=-999; phiG2_gen=-999; eG2_gen=-999;
    ptB1_gen=-999; etaB1_gen=-999; phiB1_gen=-999; mB1_gen=-999;
    ptB2_gen=-999; etaB2_gen=-999; phiB2_gen=-999; mB2_gen=-999;
    ptB3_gen=-999; etaB3_gen=-999; phiB3_gen=-999; mB3_gen=-999;
    ptB4_gen=-999; etaB4_gen=-999; phiB4_gen=-999; mB4_gen=-999;

    ptTau1_genJet=-999; etaTau1_genJet=-999; phiTau1_genJet=-999; mTau1_genJet=-999;
    ptTau2_genJet=-999; etaTau2_genJet=-999; phiTau2_genJet=-999; mTau2_genJet=-999;

    eventType=-1;

    dEta_tt=-999; dEta_6j=-999; 
    nBtag=0; nCentral=0; nJets=0; centB=0; nLep=0;
    iHmatch1=0; iHmatch2=0; iHmatch3=0; iHmatch4=0;

    isBBTT=0; isBBGG=0; isBBBB=0;
    isVBFTT=0; isVBF4B=0;

    jetTau1=0; jetTau2=0; eleTau=0;  muTau=0; 
    isoTau1=0; isoTau2=0;
    jetB1=0;   jetB2=0;   jetB3=0;   jetB4=0;
    jet_tt1=0; jet_tt2=0; jet_6j1=0; jet_6j2=0;
    gamma1=0;  gamma2=0;

    met=0; metPhi=0; ppMet=0; ppMetPhi=0; pileupMet=0; pileupMetPhi=0;
    m_sv = -999; m_svpileup = -999; m_svpuppi = -999;
    mt2=-999; mt2puppi=-999; mt2pileup=-999;
    
    rho_0=-999; rho_1=-999; rho_2=-999;

    dEtaBB1=-999; dPhiBB1=-999; dRBB1=-999;
    dEtaBB2=-999; dPhiBB2=-999; dRBB2=-999;
    dEtaB1B2=-999; dPhiB1B2=-999; dRB1B2=-999;
    dEtaB1B3=-999; dPhiB1B3=-999; dRB1B3=-999;
    dEtaB1B4=-999; dPhiB1B4=-999; dRB1B4=-999;
    dEtaB2B3=-999; dPhiB2B3=-999; dRB2B3=-999;
    dEtaB2B4=-999; dPhiB2B4=-999; dRB2B4=-999;
    dEtaB3B4=-999; dPhiB3B4=-999; dRB3B4=-999;
    dEtaTT=-999; dPhiTT=-999; dRTT=-999;
    dEtaHH=-999; dPhiHH=-999; dRHH=-999;

    mindR4B=-999;
    nBJetsComb=0;

    // ********************
    // EVENT WEIGHT
    // ********************

    eventWeight = 1;
    if (branchEvent && sampleNo!=100) {
      event = (LHEFEvent*) branchEvent->At(0);
      eventWeight*=event->Weight;
    }
    eventWeight *= xsec;

    // **********************
    // Event rho
    // *********************
    
    rho0 = (Rho*) branchRho->At(0);
    rho1 = (Rho*) branchRho->At(1);
    rho2 = (Rho*) branchRho->At(2);

    rho_0 = rho0->Rho;
    rho_1 = rho1->Rho;
    rho_2 = rho2->Rho;

    // ********************
    // TAU SELECTION
    // ********************
    for (Int_t iJet=0; iJet<branchJet->GetEntries(); iJet++) { // reconstructed jet loop
      jet = (Jet*) branchJet->At(iJet);   
      
      double corr = 1.0; //doJetcorr(corrector,jet,rho_2,rho_1,rho_0);
      
      if (fabs(jet->Eta)>4.0) continue;
      if (corr*jet->PT<20) continue;
      if (jet->TauTag==0) continue;
      
      if (puJetID(jet->Eta, jet->MeanSqDeltaR, jet->BetaStar)==1) continue;

      Int_t lep=0;
      for (Int_t iMuon=0; iMuon<branchMuon->GetEntries(); iMuon++) { // reco muon loop
	mu = (Muon*) branchMuon->At(iMuon);
	if (mu->PT>20 && deltaR(mu->Eta, jet->Eta, mu->Phi, jet->Phi)<MAX_MATCH_DIST)lep++;
      }
      for (Int_t iEle=0; iEle<branchElectron->GetEntries(); iEle++) { // reco ele loop
	ele = (Electron*) branchElectron->At(iEle);
	if (ele->PT>20 && deltaR(ele->Eta, jet->Eta, ele->Phi, jet->Phi)<MAX_MATCH_DIST)lep++;
      }
      if (lep>0) continue;
      
      if (iT1==-1) { 
	iT1=iJet; 
	jetTau1 = (Jet*) branchJet->At(iT1); 
	tauCat1=hadron;
      }
      else if (corr*jet->PT > jetTau1->PT) { 
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
      else if (corr*jet->PT > jetTau2->PT) { 
	iT2=iJet; 
	jetTau2 = (Jet*) branchJet->At(iT2); 
	tauCat2=hadron;
      }
    } // end reco jet loop

    for (Int_t iMuon=0; iMuon<branchMuon->GetEntries(); iMuon++) { // reco muon loop
      mu = (Muon*) branchMuon->At(iMuon);
      
      if (fabs(mu->Eta)>4.0) continue;
      if (mu->PT<20) continue;
      if (mu->IsolationVar>0.4) continue;

      nLep++;

      if ((jetTau1)&&(deltaR(mu->Eta, jetTau1->Eta, mu->Phi, jetTau1->Phi) < MAX_MATCH_DIST)) continue;
      if ((jetTau2)&&(deltaR(mu->Eta, jetTau2->Eta, mu->Phi, jetTau2->Phi) < MAX_MATCH_DIST)) continue;      

      if(!muTau)
	muTau = (Muon*) branchMuon->At(iMuon);
      else
	{
	  if ( mu->PT > muTau->PT ) { 
	    muTau = (Muon*) branchMuon->At(iMuon); 
	  }
	}
    }
    
    // get electronic taus
    for (Int_t iEle=0; iEle<branchElectron->GetEntries(); iEle++) { // reco ele loop
      ele = (Electron*) branchElectron->At(iEle);
      
      if (fabs(ele->Eta)>4.0) continue;
      if (ele->PT<20) continue;
      if (ele->IsolationVar>0.4) continue;

      nLep++;

      if ((jetTau1)&&(deltaR(ele->Eta, jetTau1->Eta, ele->Phi, jetTau1->Phi) < MAX_MATCH_DIST)) continue;
      if ((jetTau2)&&(deltaR(ele->Eta, jetTau2->Eta, ele->Phi, jetTau2->Phi) < MAX_MATCH_DIST)) continue;      

      if(!eleTau)
	eleTau= (Electron*) branchElectron->At(iEle);
      else
	{
	  if ( ele->PT > eleTau->PT ) { 
	    eleTau = (Electron*) branchElectron->At(iEle); 
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
      if (gam->IsolationVar>0.4) continue;

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

    //cout << "jet branch has " << branchJet->GetEntries() << endl;

    for (Int_t iJet=0; iJet<branchJet->GetEntries(); iJet++) { // reconstructed jet loop
      jet = (Jet*) branchJet->At(iJet);

      //cout << "jet " << iJet << endl;

      double corr = 1.0; //doJetcorr(corrector,jet,rho_2,rho_1,rho_0);

      if (fabs(jet->Eta)>4.0) continue;
      if (corr*jet->PT<20) continue;

      if (puJetID(jet->Eta, jet->MeanSqDeltaR, jet->BetaStar)==1) continue;
      
      nJets++;

      if (jet->BTag==0) continue;
      
      if ((jetTau1)&&(deltaR(jet->Eta, jetTau1->Eta, jet->Phi, jetTau1->Phi) < MAX_MATCH_DIST)) continue;
      if ((jetTau2)&&(deltaR(jet->Eta, jetTau2->Eta, jet->Phi, jetTau2->Phi) < MAX_MATCH_DIST)) continue;
      if ((muTau)&&(deltaR(jet->Eta, muTau->Eta, jet->Phi, muTau->Phi) < MAX_MATCH_DIST)) continue;
      if ((eleTau)&&(deltaR(jet->Eta, eleTau->Eta, jet->Phi, eleTau->Phi) < MAX_MATCH_DIST)) continue;

      if (iB1==-1) {
	iB1=iJet; 
	jetB1 = (Jet*) branchJet->At(iB1); 
	bTag1=jetB1->BTag;
      }
      else if (corr*jet->PT > jetB1->PT) {
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
      else if (corr*jet->PT > jetB2->PT) { 
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
      else if (corr*jet->PT > jetB3->PT) {
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
      else if (corr*jet->PT > jetB4->PT) {
	iB4=iJet;
	jetB4 = (Jet*) branchJet->At(iB4);
	bTag4=jetB4->BTag;
      }
    }

    // ********************
    // VBF-JET SELECTION
    // ********************
    for (Int_t iJet=0; iJet<branchRawJet->GetEntries(); iJet++) { // reconstructed jet loop
      jet = (Jet*) branchRawJet->At(iJet);
      double corr = doJetcorr(corrector,jet,rho_2,rho_1,rho_0);
       
      if (fabs(jet->Eta)>4.7) continue;
      if (corr*jet->PT<30) continue;

      if (puJetID(jet->Eta, jet->MeanSqDeltaR, jet->BetaStar)==1) continue;

      if (jet->BTag>0) nBtag++;

      if ((jetTau1)&&(deltaR(jet->Eta, jetTau1->Eta, jet->Phi, jetTau1->Phi) < MAX_MATCH_DIST)) continue;
      if ((jetTau2)&&(deltaR(jet->Eta, jetTau2->Eta, jet->Phi, jetTau2->Phi) < MAX_MATCH_DIST)) continue;
      if ((muTau)&&(deltaR(jet->Eta, muTau->Eta, jet->Phi, muTau->Phi) < MAX_MATCH_DIST)) continue;
      if ((eleTau)&&(deltaR(jet->Eta, eleTau->Eta, jet->Phi, eleTau->Phi) < MAX_MATCH_DIST)) continue;

      if (iJ1==-1) {
        iJ1=iJet;
        jet_tt1 = (Jet*) branchRawJet->At(iJ1);
        jbTag_tt1=jet_tt1->BTag;
      }
      else if (corr*jet->PT >  doJetcorr(corrector,jet_tt1,rho_2,rho_1,rho_0)*jet_tt1->PT) {
        iJ2=iJ1;
        jet_tt2 = (Jet*) branchRawJet->At(iJ2);
        jbTag_tt2=jet_tt2->BTag;
        iJ1=iJet;
        jet_tt1 = (Jet*) branchRawJet->At(iJ1);
        jbTag_tt1=jet_tt1->BTag;
      }
      else if (iJ2==-1) {
        iJ2=iJet;
        jet_tt2 = (Jet*) branchRawJet->At(iJ2);
        jbTag_tt2=jet_tt2->BTag;
      }
      else if (corr*jet->PT >  doJetcorr(corrector,jet_tt2,rho_2,rho_1,rho_0)*jet_tt2->PT) {
        iJ2=iJet;
        jet_tt2 = (Jet*) branchRawJet->At(iJ2);
        jbTag_tt2=jet_tt2->BTag;
      }
    }
    for (Int_t iJet=0; iJet<branchRawJet->GetEntries(); iJet++) { // reconstructed jet loop
      jet = (Jet*) branchRawJet->At(iJet);
      double corr = doJetcorr(corrector,jet,rho_2,rho_1,rho_0);
      
      if (fabs(jet->Eta)>4.7) continue;
      if (corr*jet->PT<30) continue;
      if(fabs(jet->Eta)>4.0) 

      if (puJetID(jet->Eta, jet->MeanSqDeltaR, jet->BetaStar)==1) continue;

      if ((jetB1)&&(deltaR(jet->Eta, jetB1->Eta, jet->Phi, jetB1->Phi) < MAX_MATCH_DIST)) continue;
      if ((jetB2)&&(deltaR(jet->Eta, jetB2->Eta, jet->Phi, jetB2->Phi) < MAX_MATCH_DIST)) continue;
      if ((jetB3)&&(deltaR(jet->Eta, jetB3->Eta, jet->Phi, jetB3->Phi) < MAX_MATCH_DIST)) continue;
      if ((jetB4)&&(deltaR(jet->Eta, jetB4->Eta, jet->Phi, jetB4->Phi) < MAX_MATCH_DIST)) continue;

      if (iJ3==-1) {
        iJ3=iJet;
        jet_6j1 = (Jet*) branchRawJet->At(iJ3);
        jbTag_6j1=jet_6j1->BTag;
      }
      else if (corr*jet->PT > doJetcorr(corrector,jet_6j1,rho_2,rho_1,rho_0)*jet_6j1->PT) {
	iJ4=iJ3;
	jet_6j2 = (Jet*) branchRawJet->At(iJ4);
        jbTag_6j2=jet_6j2->BTag;
        iJ3=iJet;
        jet_6j1 = (Jet*) branchRawJet->At(iJ3);
        jbTag_6j1=jet_6j1->BTag;
      }
      else if (iJ4==-1) {
        iJ4=iJet;
        jet_6j2 = (Jet*) branchRawJet->At(iJ4);
        jbTag_6j2=jet_6j2->BTag;
      }
      else if (corr*jet->PT > doJetcorr(corrector,jet_6j2,rho_2,rho_1,rho_0)*jet_6j2->PT) {
        iJ4=iJet;
        jet_6j2 = (Jet*) branchRawJet->At(iJ4);
        jbTag_6j2=jet_6j2->BTag;
      }
    }

    // ********************
    // STORE VARIABLES
    // ********************

    //std::cout << "whaaat  " <<  (muTau ? muTau->PT : 0) << "  "  << (jetTau1 ? jetTau1->PT : 0) << " "  << (jetTau2 ? jetTau2->PT : 0) << "  " << (eleTau ? eleTau->PT : 0) << std::endl;
    // fill 4-vector for leading tau
    LorentzVector vRecoTau1(0,0,0,0);
    LorentzVector vRecoTau2(0,0,0,0);
    double corrtau1=0, corrtau2=0;
    if(jetTau1)
      corrtau1 = 1.0; //doJetcorr(corrector,jetTau1,rho_2,rho_1,rho_0);
    if(jetTau2)
      corrtau2 = 1.0; //doJetcorr(corrector,jetTau2,rho_2,rho_1,rho_0);
    if(jetTau2 ? (corrtau2*jetTau2->PT > (muTau ? muTau->PT : 0) && corrtau2*jetTau2->PT > (eleTau? eleTau->PT : 0)) : 0) {
      vRecoTau1.SetPt(corrtau1*jetTau1->PT);
      vRecoTau1.SetEta(jetTau1->Eta);
      vRecoTau1.SetPhi(jetTau1->Phi);
      vRecoTau1.SetM(corrtau1*jetTau1->Mass);
      ptTau1=corrtau1*jetTau1->PT;
      etaTau1=jetTau1->Eta;
      phiTau1=jetTau1->Phi;
      mTau1=corrtau1*jetTau1->Mass;
      vRecoTau2.SetPt(corrtau2*jetTau2->PT);
      vRecoTau2.SetEta(jetTau2->Eta);
      vRecoTau2.SetPhi(jetTau2->Phi);
      vRecoTau2.SetM(corrtau2*jetTau2->Mass);
      ptTau2=corrtau2*jetTau2->PT;
      etaTau2=jetTau2->Eta;
      phiTau2=jetTau2->Phi;
      mTau2=corrtau2*jetTau2->Mass;
    }
    else if ((muTau ? muTau->PT : 0) > (jetTau1 ? jetTau1->PT : 0) && (eleTau ? eleTau->PT : 0) > (jetTau1 ? jetTau1->PT : 0)) {
      vRecoTau1.SetPt(muTau->PT);
      vRecoTau1.SetEta(muTau->Eta);
      vRecoTau1.SetPhi(muTau->Phi);
      vRecoTau1.SetM(MUON_MASS);
      ptTau1=muTau->PT;
      etaTau1=muTau->Eta;
      phiTau1=muTau->Phi;
      mTau1=MUON_MASS;
      tauCat1=muon;
      tauIso1=muTau->IsolationVar;
      vRecoTau2.SetPt(eleTau->PT);
      vRecoTau2.SetEta(eleTau->Eta);
      vRecoTau2.SetPhi(eleTau->Phi);
      vRecoTau2.SetM(ELE_MASS);
      ptTau2=eleTau->PT;
      etaTau2=eleTau->Eta;
      phiTau2=eleTau->Phi;
      mTau2=ELE_MASS;
      tauCat2=electron;
      tauIso2=eleTau->IsolationVar;
    }
    else if (jetTau1 && (muTau? muTau->PT : 0) > (jetTau2 ? corrtau2*jetTau2->PT : 0) && (muTau ? muTau->PT : 0) > (eleTau ? eleTau->PT : 0)) {
      vRecoTau1.SetPt(corrtau1*jetTau1->PT);
      vRecoTau1.SetEta(jetTau1->Eta);
      vRecoTau1.SetPhi(jetTau1->Phi);
      vRecoTau1.SetM(corrtau1*jetTau1->Mass);
      ptTau1=corrtau1*jetTau1->PT;
      etaTau1=jetTau1->Eta;
      phiTau1=jetTau1->Phi;
      mTau1=corrtau1*jetTau1->Mass;
      vRecoTau2.SetPt(muTau->PT);
      vRecoTau2.SetEta(muTau->Eta);
      vRecoTau2.SetPhi(muTau->Phi);
      vRecoTau2.SetM(MUON_MASS);
      ptTau2=muTau->PT;
      etaTau2=muTau->Eta;
      phiTau2=muTau->Phi;
      mTau2=MUON_MASS;
      tauCat2=muon;
      tauIso2=muTau->IsolationVar;
    }
    else if (jetTau1 && (eleTau ? eleTau->PT : 0) > (jetTau2 ? corrtau2*jetTau2->PT : 0) && (eleTau ? eleTau->PT : 0) > (muTau ? muTau->PT : 0)) {
      vRecoTau1.SetPt(corrtau1*jetTau1->PT);
      vRecoTau1.SetEta(jetTau1->Eta);
      vRecoTau1.SetPhi(jetTau1->Phi);
      vRecoTau1.SetM(corrtau1*jetTau1->Mass);
      ptTau1=corrtau1*jetTau1->PT;
      etaTau1=jetTau1->Eta;
      phiTau1=jetTau1->Phi;
      mTau1=corrtau1*jetTau1->Mass;
      vRecoTau2.SetPt(eleTau->PT);
      vRecoTau2.SetEta(eleTau->Eta);
      vRecoTau2.SetPhi(eleTau->Phi);
      vRecoTau2.SetM(ELE_MASS);
      ptTau2=eleTau->PT;
      etaTau2=eleTau->Eta;
      phiTau2=eleTau->Phi;
      mTau2=ELE_MASS;
      tauCat2=electron;
      tauIso2=eleTau->IsolationVar;
    }
    
    if (tauCat1==hadron) {
      for (Int_t iIso=0; iIso<branchIsoTrack->GetEntries(); iIso++) {
	iso=(IsoTrack*) branchIsoTrack->At(iIso);

	if (iso->IsolationVar>0.4) continue;

	if (deltaR(iso->Eta, etaTau1, iso->Phi, phiTau1)<MAX_MATCH_DIST) {
	  if (!(isoTau1)) {
	    iI1=iIso;
	    isoTau1=(IsoTrack*) branchIsoTrack->At(iI1);
	  }
	  else if (isoTau1 && iso->PT>isoTau1->PT) {
	    iI1=iIso;
	    isoTau1=(IsoTrack*) branchIsoTrack->At(iI1);
	  }
	}
      }
    }

    if (tauCat2==hadron) {
      for (Int_t iIso=0; iIso<branchIsoTrack->GetEntries(); iIso++) {
	iso=(IsoTrack*) branchIsoTrack->At(iIso);
	
	if (iso->IsolationVar>0.4) continue;

	if (deltaR(iso->Eta, etaTau2, iso->Phi, phiTau2)<MAX_MATCH_DIST) {
	  if (!(isoTau2)) {
	    iI2=iIso;
	    isoTau2=(IsoTrack*) branchIsoTrack->At(iI2);
	  }
	  else if (isoTau2 && iso->PT>isoTau2->PT) {
	    iI2=iIso;
	    isoTau2=(IsoTrack*) branchIsoTrack->At(iI2);
	  }
	}
      }
    }

    if (isoTau1) {
      ptTrk1=isoTau1->PT;
      etaTrk1=isoTau1->Eta;
      phiTrk1=isoTau1->Phi;
      mTrk1=0;
    }
    if (isoTau2) {
      ptTrk2=isoTau2->PT;
      etaTrk2=isoTau2->Eta;
      phiTrk2=isoTau2->Phi;
      mTrk2=0;
    }

    // fill 4-vector for leading b-jet
    LorentzVector vRecoB1(0,0,0,0);
    if (jetB1) {
      double corrb1 = 1.0; //doJetcorr(corrector,jetB1,rho_2,rho_1,rho_0);
      vRecoB1.SetPt(corrb1*jetB1->PT);
      vRecoB1.SetEta(jetB1->Eta);
      vRecoB1.SetPhi(jetB1->Phi);
      vRecoB1.SetM(corrb1*jetB1->Mass);
      ptB1=corrb1*jetB1->PT;
      etaB1=jetB1->Eta;
      phiB1=jetB1->Phi;
      mB1=corrb1*jetB1->Mass;
    }

    // fill 4-vector for second b-jet
    LorentzVector vRecoB2(0,0,0,0);
    if (jetB2) {
      double corrb2 = 1.0; //doJetcorr(corrector,jetB2,rho_2,rho_1,rho_0);
      vRecoB2.SetPt(corrb2*jetB2->PT);
      vRecoB2.SetEta(jetB2->Eta);
      vRecoB2.SetPhi(jetB2->Phi);
      vRecoB2.SetM(corrb2*jetB2->Mass);
      ptB2=corrb2*jetB2->PT;
      etaB2=jetB2->Eta;
      phiB2=jetB2->Phi;
      mB2=corrb2*jetB2->Mass;
    }

    // fill 4-vector for b-jet
    LorentzVector vRecoB3(0,0,0,0);
    if (jetB3) {
      double corrb3 = 1.0; //doJetcorr(corrector,jetB3,rho_2,rho_1,rho_0);
      vRecoB3.SetPt(corrb3*jetB3->PT);
      vRecoB3.SetEta(jetB3->Eta);
      vRecoB3.SetPhi(jetB3->Phi);
      vRecoB3.SetM(corrb3*jetB3->Mass);
      ptB3=corrb3*jetB3->PT;
      etaB3=jetB3->Eta;
      phiB3=jetB3->Phi;
      mB3=corrb3*jetB3->Mass;
    }
    // fill 4-vector for b-jet
    LorentzVector vRecoB4(0,0,0,0);
    if (jetB4) {
      double corrb4 = 1.0; //doJetcorr(corrector,jetB4,rho_2,rho_1,rho_0);
      vRecoB4.SetPt(corrb4*jetB4->PT);
      vRecoB4.SetEta(jetB4->Eta);
      vRecoB4.SetPhi(jetB4->Phi);
      vRecoB4.SetM(corrb4*jetB4->Mass);
      ptB4=corrb4*jetB4->PT;
      etaB4=jetB4->Eta;
      phiB4=jetB4->Phi;
      mB4=corrb4*jetB4->Mass;
    }

    // fill 4-vector for leading photon
    LorentzVector vRecoGam1(0,0,0,0);
    if (gamma1) {
      vRecoGam1.SetPt(gamma1->PT);
      vRecoGam1.SetEta(gamma1->Eta);
      vRecoGam1.SetPhi(gamma1->Phi);
      vRecoGam1.SetM(0);
      ptG1=gamma1->PT;
      etaG1=gamma1->Eta;
      phiG1=gamma1->Phi;
      eG1=gamma1->E;
    }

    // fill 4-vector for second photon
    LorentzVector vRecoGam2(0,0,0,0);
    if (gamma2) {
      vRecoGam2.SetPt(gamma2->PT);
      vRecoGam2.SetEta(gamma2->Eta);
      vRecoGam2.SetPhi(gamma2->Phi);
      vRecoGam2.SetM(0);
      ptG2=gamma2->PT;
      etaG2=gamma2->Eta;
      phiG2=gamma2->Phi;
      eG2=gamma2->E;
    }
    // fill 4-vector for leading VBF jet
    LorentzVector vRecoJet_tt1(0,0,0,0);
    if (jet_tt1) {
      double corrj1 = doJetcorr(corrector,jet_tt1,rho_2,rho_1,rho_0);
      vRecoJet_tt1.SetPt(corrj1*jet_tt1->PT);
      vRecoJet_tt1.SetEta(jet_tt1->Eta);
      vRecoJet_tt1.SetPhi(jet_tt1->Phi);
      vRecoJet_tt1.SetM(corrj1*jet_tt1->Mass);
      ptJet_tt1=corrj1*jet_tt1->PT;
      etaJet_tt1=jet_tt1->Eta;
      phiJet_tt1=jet_tt1->Phi;
      mJet_tt1=corrj1*jet_tt1->Mass;
    }

    LorentzVector vRecoJet_tt2(0,0,0,0);
    if (jet_tt2) {
      double corrj2 = doJetcorr(corrector,jet_tt2,rho_2,rho_1,rho_0);
      vRecoJet_tt2.SetPt(corrj2*jet_tt2->PT);
      vRecoJet_tt2.SetEta(jet_tt2->Eta);
      vRecoJet_tt2.SetPhi(jet_tt2->Phi);
      vRecoJet_tt2.SetM(corrj2*jet_tt2->Mass);
      ptJet_tt2=corrj2*jet_tt2->PT;
      etaJet_tt2=jet_tt2->Eta;
      phiJet_tt2=jet_tt2->Phi;
      mJet_tt2=corrj2*jet_tt2->Mass;
    }
    // fill 4-vector for leading VBF jet
    LorentzVector vRecoJet_6j1(0,0,0,0);
    if (jet_6j1) {
      double corrj1 = doJetcorr(corrector,jet_6j1,rho_2,rho_1,rho_0);
      vRecoJet_6j1.SetPt(corrj1*jet_6j1->PT);
      vRecoJet_6j1.SetEta(jet_6j1->Eta);
      vRecoJet_6j1.SetPhi(jet_6j1->Phi);
      vRecoJet_6j1.SetM(corrj1*jet_6j1->Mass);
      ptJet_6j1=corrj1*jet_6j1->PT;
      etaJet_6j1=jet_6j1->Eta;
      phiJet_6j1=jet_6j1->Phi;
      mJet_6j1=corrj1*jet_6j1->Mass;
    }

    LorentzVector vRecoJet_6j2(0,0,0,0);
    if (jet_6j2) {
      double corrj2 = doJetcorr(corrector,jet_6j2,rho_2,rho_1,rho_0);
      vRecoJet_6j2.SetPt(corrj2*jet_6j2->PT);
      vRecoJet_6j2.SetEta(jet_6j2->Eta);
      vRecoJet_6j2.SetPhi(jet_6j2->Phi);
      vRecoJet_6j2.SetM(corrj2*jet_6j2->Mass);
      ptJet_6j2=corrj2*jet_6j2->PT;
      etaJet_6j2=jet_6j2->Eta;
      phiJet_6j2=jet_6j2->Phi;
      mJet_6j2=corrj2*jet_6j2->Mass;
    }
    // ********************
    // FLAG FINAL STATE
    // ********************

    if (ptB1>0 && ptB2>0 && ptTau1>0 && ptTau2>0) { isBBTT=1; }
    if (ptB1>0 && ptB2>0 && ptB3>0 && ptB4>0) { isBBBB=1; }
    if (ptB1>0 && ptB2>0 && ptG1>0 && ptG2>0) { isBBGG=1; }
    if (ptJet_tt1>0 && ptJet_tt2>0 && ptTau1>0 && ptTau2>0) { isVBFTT=1; }
    if (ptB1>0 && ptB2>0 && ptB3>0 && ptB4>0 && ptJet_6j1>0 && ptJet_6j2>0) { isVBF4B=1; }
    if (isBBTT==0 && isBBBB==0 && isBBGG==0 && isVBFTT==0 && isVBF4B==0) continue;

    
    // ********************
    // COMPUTE VARIABLES
    // ********************

    LorentzVector vTT;
    if (isBBTT==1 || isVBFTT==1) {
      vTT = vRecoTau1+vRecoTau2;
      mTT=vTT.M();
      ptTT=vTT.Pt();
      etaTT=vTT.Eta();
      phiTT=vTT.Phi();
      dEtaTT=abs(vRecoTau1.Eta()-vRecoTau2.Eta());
      dPhiTT=deltaPhi(vRecoTau1.Phi(), vRecoTau2.Phi());
      dRTT=deltaR(vRecoTau1.Eta(), vRecoTau2.Eta(), vRecoTau1.Phi(), vRecoTau2.Phi());
    }

    LorentzVector vGG;
    if (isBBGG==1) {
      vGG = vRecoGam1+vRecoGam2;
      mGG=vGG.M();
      ptGG=vGG.Pt();
      etaGG=vGG.Eta();
      phiGG=vGG.Phi();
    }
    
    LorentzVector vBB1;
    if (isBBTT==1 || isBBGG==1 || isBBBB==1 || isVBF4B==1) {
      vBB1 = vRecoB1+vRecoB2;
      mBB1=vBB1.M();
      ptBB1=vBB1.Pt();
      etaBB1=vBB1.Eta();
      phiBB1=vBB1.Phi();
      dEtaBB1=abs(vRecoB1.Eta()-vRecoB2.Eta());
      dPhiBB1=deltaPhi(vRecoB1.Phi(), vRecoB2.Phi());
      dRBB1=deltaR(vRecoB1.Eta(), vRecoB2.Eta(), vRecoB1.Phi(), vRecoB2.Phi());
    }

    LorentzVector vBB2;
    if (isBBBB==1 || isVBF4B==1) {
      vBB2 = vRecoB3+vRecoB4;
      mBB2=vBB2.M();
      ptBB2=vBB2.Pt();
      etaBB2=vBB2.Eta();
      phiBB2=vBB2.Phi();
      dEtaBB2=abs(vRecoB3.Eta()-vRecoB4.Eta());
      dPhiBB2=deltaPhi(vRecoB3.Phi(), vRecoB4.Phi());
      dRBB2=deltaR(vRecoB3.Eta(), vRecoB4.Eta(), vRecoB3.Phi(), vRecoB4.Phi());
    }

    LorentzVector vB1B2;
    if (isBBBB==1 || isVBF4B==1) {
      vB1B2 = vRecoB1+vRecoB2;
      mB1B2=vB1B2.M();
      ptB1B2=vB1B2.Pt();
      etaB1B2=vB1B2.Eta();
      phiB1B2=vB1B2.Phi();
      dEtaB1B2=abs(vRecoB1.Eta()-vRecoB2.Eta());
      dPhiB1B2=deltaPhi(vRecoB1.Phi(), vRecoB2.Phi());
      dRB1B2=deltaR(vRecoB1.Eta(), vRecoB2.Eta(), vRecoB1.Phi(), vRecoB2.Phi());
    }
    
    LorentzVector vB1B3;
    if (isBBBB==1 || isVBF4B==1) {
      vB1B3 = vRecoB1+vRecoB3;
      mB1B3=vB1B3.M();
      ptB1B3=vB1B3.Pt();
      etaB1B3=vB1B3.Eta();
      phiB1B3=vB1B3.Phi();
      dEtaB1B3=abs(vRecoB1.Eta()-vRecoB3.Eta());
      dPhiB1B3=deltaPhi(vRecoB1.Phi(), vRecoB3.Phi());
      dRB1B3=deltaR(vRecoB1.Eta(), vRecoB3.Eta(), vRecoB1.Phi(), vRecoB3.Phi());
    }

    LorentzVector vB1B4;
    if (isBBBB==1 || isVBF4B==1) {
      vB1B4 = vRecoB1+vRecoB4;
      mB1B4=vB1B4.M();
      ptB1B4=vB1B4.Pt();
      etaB1B4=vB1B4.Eta();
      phiB1B4=vB1B4.Phi();
      dEtaB1B4=abs(vRecoB1.Eta()-vRecoB4.Eta());
      dPhiB1B4=deltaPhi(vRecoB1.Phi(), vRecoB4.Phi());
      dRB1B4=deltaR(vRecoB1.Eta(), vRecoB4.Eta(), vRecoB1.Phi(), vRecoB4.Phi());
    }

    LorentzVector vB2B3;
    if (isBBBB==1 || isVBF4B==1) {
      vB2B3 = vRecoB2+vRecoB3;
      mB2B3=vB2B3.M();
      ptB2B3=vB2B3.Pt();
      etaB2B3=vB2B3.Eta();
      phiB2B3=vB2B3.Phi();
      dEtaB2B3=abs(vRecoB2.Eta()-vRecoB3.Eta());
      dPhiB2B3=deltaPhi(vRecoB2.Phi(), vRecoB3.Phi());
      dRB2B3=deltaR(vRecoB2.Eta(), vRecoB3.Eta(), vRecoB2.Phi(), vRecoB3.Phi());
    }

    LorentzVector vB2B4;
    if (isBBBB==1 || isVBF4B==1) {
      vB2B4 = vRecoB2+vRecoB4;
      mB2B4=vB2B4.M();
      ptB2B4=vB2B4.Pt();
      etaB2B4=vB2B4.Eta();
      phiB2B4=vB2B4.Phi();
      dEtaB2B4=abs(vRecoB2.Eta()-vRecoB4.Eta());
      dPhiB2B4=deltaPhi(vRecoB2.Phi(), vRecoB4.Phi());
      dRB2B4=deltaR(vRecoB2.Eta(), vRecoB4.Eta(), vRecoB2.Phi(), vRecoB4.Phi());
    }

    LorentzVector vB3B4;
    if (isBBBB==1 || isVBF4B==1) {
      vB3B4 = vRecoB3+vRecoB4;
      mB3B4=vB3B4.M();
      ptB3B4=vB3B4.Pt();
      etaB3B4=vB3B4.Eta();
      phiB3B4=vB3B4.Phi();
      dEtaB3B4=abs(vRecoB3.Eta()-vRecoB4.Eta());
      dPhiB3B4=deltaPhi(vRecoB3.Phi(), vRecoB4.Phi());
      dRB3B4=deltaR(vRecoB3.Eta(), vRecoB4.Eta(), vRecoB3.Phi(), vRecoB4.Phi());
    }

    if (isBBBB==1 || isVBF4B==1) {
      mindR4B=min(min(min(min(min(dRB1B2,dRB1B3),dRB1B4),dRB2B3),dRB2B4),dRB3B4);
      if(90<mB1B2<135)
	{
	  nBJetsComb++;
	}
      if(90<mB1B3<135)
	{
	  nBJetsComb++;
	}
      if(90<mB1B4<135)
	{
	  nBJetsComb++;
	}
      if(90<mB2B3<135)
	{
	  nBJetsComb++;
	}
      if(90<mB2B4<135)
	{
	  nBJetsComb++;
	}
      if(90<mB3B4<135)
	{
	  nBJetsComb++;
	}
    }

    LorentzVector vHH;
    if (isBBTT==1) {
      vHH = vTT+vBB1;
      mHH=vHH.M();
      ptHH=vHH.Pt();
      etaHH=vHH.Eta();
      phiHH=vHH.Phi();
      dEtaHH=abs(vBB1.Eta()-vTT.Eta());
      dPhiHH=deltaPhi(vBB1.Phi(), vTT.Phi());
      dRHH=deltaR(vBB1.Eta(), vTT.Eta(), vBB1.Phi(), vTT.Phi());
    }
    else if (isBBGG==1) {
      vHH = vTT+vGG;
      mHH=vHH.M();
      ptHH=vHH.Pt();
      etaHH=vHH.Eta();
      phiHH=vHH.Phi();
    }
    else if (isBBBB==1 || isVBF4B==1) {
      vHH = vBB1+vBB2;
      mHH=vHH.M();
      ptHH=vHH.Pt();
      etaHH=vHH.Eta();
      phiHH=vHH.Phi();
    }
    
    LorentzVector vJJ;
    if (vRecoJet_tt1.Pt()>0 && vRecoJet_tt2.Pt()>0) {
      vJJ = vRecoJet_tt1+vRecoJet_tt2;
      mJJ_tt=vJJ.M();
      ptJJ_tt=vJJ.Pt();
      etaJJ_tt=vJJ.Eta();
      phiJJ_tt=vJJ.Phi();
      dEta_tt=vRecoJet_tt1.Eta()-vRecoJet_tt2.Eta();
     
      for (Int_t iJet=0; iJet<branchRawJet->GetEntries(); iJet++) { // reconstructed jet loop                                                                                                  
	jet = (Jet*) branchRawJet->At(iJet);
	double corr = doJetcorr(corrector,jet,rho_2,rho_1,rho_0);

	if (fabs(jet->Eta)>4.7) continue;
	if (corr*jet->PT<30) continue;
	if (puJetID(jet->Eta, jet->MeanSqDeltaR, jet->BetaStar)==1) continue;
	
	if ((jetTau1)&&(deltaR(jet->Eta, jetTau1->Eta, jet->Phi, jetTau1->Phi) < MAX_MATCH_DIST)) continue;
	if ((jetTau2)&&(deltaR(jet->Eta, jetTau2->Eta, jet->Phi, jetTau2->Phi) < MAX_MATCH_DIST)) continue;
	if ((muTau)&&(deltaR(jet->Eta, muTau->Eta, jet->Phi, muTau->Phi) < MAX_MATCH_DIST)) continue;
	if ((eleTau)&&(deltaR(jet->Eta, eleTau->Eta, jet->Phi, eleTau->Phi) < MAX_MATCH_DIST)) continue;
	
	if ( (vRecoJet_tt1.Eta() > vRecoJet_tt2.Eta()) && (jet->Eta > vRecoJet_tt2.Eta()) && (vRecoJet_tt1.Eta() > jet->Eta) ) {
	  nCentral++;
	}
	else if ( (vRecoJet_tt2.Eta() > vRecoJet_tt1.Eta()) && (jet->Eta > vRecoJet_tt1.Eta()) && (vRecoJet_tt2.Eta() > jet->Eta) ) {
	  nCentral++;
	}
      }	
    }

    if (vRecoJet_6j1.Pt()>0 && vRecoJet_6j2.Pt()>0) {
      vJJ = vRecoJet_6j1+vRecoJet_6j2;
      mJJ_6j=vJJ.M();
      ptJJ_6j=vJJ.Pt();
      etaJJ_6j=vJJ.Eta();
      phiJJ_6j=vJJ.Phi();
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
    // ********************
    // MT2 CALC
    // ********************

    missET = (MissingET*) branchMET->At(0);

    met=missET->MET;
    metPhi=missET->Phi;

    if (0) {
    //if (branchPuppiMET) {
      missET = (MissingET*) branchPuppiMET->At(0);
      ppMet=missET->MET;
      ppMetPhi=missET->Phi;
    }

    if (branchPileupMET) {
      missET = (MissingET*) branchPileupMET->At(0);
      pileupMet=missET->MET;
      pileupMetPhi=missET->Phi;
    }

    //if (0) {
    if ( vRecoTau1.Pt()>0 && vRecoTau2.Pt()>0 && vRecoB1.Pt()>0 && vRecoB2.Pt()>0) {

      tau1.SetMagPhi(vRecoTau1.Pt(), vRecoTau1.Phi());
      tau2.SetMagPhi(vRecoTau2.Pt(), vRecoTau2.Phi());
      
      b1.SetMagPhi(vRecoB1.Pt(), vRecoB1.Phi());
      b2.SetMagPhi(vRecoB2.Pt(), vRecoB2.Phi());
      
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
      
      double step[2] = {0.1, 0.1};
      double variable[2] = { 0.5*c1.Mod(), 0.0 };
      
      ROOT::Math::Functor f(calcmt2,2);
      
      min->SetFunction(f);
      min->SetLimitedVariable(0,"cT",variable[0], step[0], 0.0, sumPt.Mod());
      min->SetLimitedVariable(1,"cPhi",variable[1], step[1], 0.0, TMath::Pi());
      
      min->Minimize();
      mt2 = min->MinValue();

      if (0) {
      //if (branchPuppiMET) {	
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
	mt2puppi = min->MinValue();
      }
      else mt2puppi = -999;
      
      if (branchPileupMET) {
	
	pumpt.SetMagPhi(pileupMet, pileupMetPhi);
	TVector2 puSumPt = tau1+tau2+pumpt;
	
	calcmt2.SetMPT(puSumPt);
	c1=puSumPt;
	c2=sumPt-c1;
	
	ROOT::Math::Functor f3(calcmt2,2);
	step[0] = 0.1; step[1] = 0.1;
	variable[0] = 0.5*c1.Mod(); variable[1] = 0.0;
	
	min->SetFunction(f3);
	min->SetLimitedVariable(0,"cT",variable[0], step[0], 0.0, puSumPt.Mod());
	min->SetLimitedVariable(1,"cPhi",variable[1], step[1], 0.0, TMath::Pi());
	
	min->Minimize();
	mt2pileup = min->MinValue();
      }
      else mt2pileup = -999;
    }    
    
    // ***********************************
    // Let's start with SVFit calculations
    // ***********************************
    
    //if (0) {
    if (vRecoTau1.Pt()>0 && vRecoTau2.Pt()>0) {
      int channel=0;
      mithep::TSVfit svfit;
      svfit.cov_00=lcov00;
      svfit.cov_01=lcov01;
      svfit.cov_10=lcov10;
      svfit.cov_11=lcov11;
      TLorentzVector lvec1;
      if(tauCat1==hadron)
	lvec1.SetPtEtaPhiM(ptTau1,etaTau1,phiTau1,TAU_MASS);
      else
	lvec1.SetPtEtaPhiM(ptTau1,etaTau1,phiTau1,mTau1);
      TLorentzVector lvec2; 
      if(tauCat2==hadron)
	lvec2.SetPtEtaPhiM(ptTau2,etaTau2,phiTau2,TAU_MASS);
      else
	lvec2.SetPtEtaPhiM(ptTau2,etaTau2,phiTau2,mTau2);
      mithep::FourVectorM svlep1; svlep1.SetPxPyPzE(lvec1.Px(),lvec1.Py(),lvec1.Pz(),lvec1.E());
      mithep::FourVectorM svlep2; svlep2.SetPxPyPzE(lvec2.Px(),lvec2.Py(),lvec2.Pz(),lvec2.E());
      if(tauCat1!=hadron)
	{
	  svfit.daughter1 = svlep1;
	  svfit.daughter2 = svlep2;
	  svfit.daughterId1 = 1;
	  svfit.daughterId2 = 2;
	  if(tauCat2!=hadron)
	    channel=0;
	  else
	    channel=1;
	}
      else
	{
	  svfit.daughter1 = svlep2;
	  svfit.daughter2 = svlep1;
	  svfit.daughterId1 = 2;
	  svfit.daughterId2 = 1;
	  if(tauCat2==hadron)
	    channel=2;
	  else
	    channel=1;
	}
      m_sv = fitter->integrate(&svfit,met,metPhi,channel);
      svfit.cov_00=lcov00pp;
      svfit.cov_01=lcov01pp;
      svfit.cov_10=lcov10pp;
      svfit.cov_11=lcov11pp;
      m_svpileup = fitter->integrate(&svfit,pileupMet,pileupMetPhi,channel);
      m_svpuppi = 0.0; // fitter->integrate(&svfit,ppMet,ppMetPhi,channel);
      std::cout << tauCat1 << "  " <<  tauCat2 << "   " << channel  << "  "  << m_sv << "  "  <<  m_svpileup <<  "  " <<  m_svpuppi << "  " << std::endl; 
    }
    
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

    if (genTau1) {
      ptTau1_gen=genTau1->PT;
      etaTau1_gen=genTau1->Eta;
      phiTau1_gen=genTau1->Phi;
      mTau1_gen=genTau1->Mass;
    }

    if (genTau2) {
      ptTau2_gen=genTau2->PT;
      etaTau2_gen=genTau2->Eta;
      phiTau2_gen=genTau2->Phi;
      mTau2_gen=genTau2->Mass;
    }

    if (genGam1) {
      ptG1_gen=genGam1->PT;
      etaG1_gen=genGam1->Eta;
      phiG1_gen=genGam1->Phi;
      eG1_gen=genGam1->E;
    }

    if (genGam2) {
      ptG2_gen=genGam2->PT;
      etaG2_gen=genGam2->Eta;
      phiG2_gen=genGam2->Phi;
      eG2_gen=genGam2->E;
    }

    if (genB1) {
      ptB1_gen=genB1->PT;
      etaB1_gen=genB1->Eta;
      phiB1_gen=genB1->Phi;
      mB1_gen=genB1->Mass;
    }

    if (genB2) {
      ptB2_gen=genB2->PT;
      etaB2_gen=genB2->Eta;
      phiB2_gen=genB2->Phi;
      mB2_gen=genB2->Mass;
    }

    if (genB3) {
      ptB3_gen=genB3->PT;
      etaB3_gen=genB3->Eta;
      phiB3_gen=genB3->Phi;
      mB3_gen=genB3->Mass;
    }

    if (genB4) {
      ptB4_gen=genB4->PT;
      etaB4_gen=genB4->Eta;
      phiB4_gen=genB4->Phi;
      mB4_gen=genB4->Mass;
    }

    if (genH1) {
      ptH1_gen=genH1->PT;
      etaH1_gen=genH1->Eta;
      phiH1_gen=genH1->Phi;
      mH1_gen=genH1->Mass;
    }

    if (genH2) {
      ptH2_gen=genH2->PT;
      etaH2_gen=genH2->Eta;
      phiH2_gen=genH2->Phi;
      mH2_gen=genH2->Mass;
    }

    // ********************
    // GEN JETS
    // ********************    

    for (Int_t iJet=0; iJet<branchGenJet->GetEntries(); iJet++) { // generator level jet loop                                                                                      
      genJet = (Jet*) branchGenJet->At(iJet);

      if (deltaR(genJet->Eta, vRecoTau1.Eta(), genJet->Phi, vRecoTau1.Phi()) < MAX_MATCH_DIST) {
        iGenJetTau1=iJet;
        genJetTau1 = (Jet*) branchGenJet->At(iGenJetTau1);
      }
      else if (deltaR(genJet->Eta, vRecoTau2.Eta(), genJet->Phi, vRecoTau2.Phi()) < MAX_MATCH_DIST) {
        iGenJetTau2=iJet;
        genJetTau2 = (Jet*) branchGenJet->At(iGenJetTau2);
      }
      else if (deltaR(genJet->Eta, vRecoJet_tt1.Eta(), genJet->Phi, vRecoJet_tt1.Phi()) < MAX_MATCH_DIST) {
        iGenJet_tt1=iJet;
        genJet_tt1 = (Jet*) branchGenJet->At(iGenJet_tt1);
      }
      else if (deltaR(genJet->Eta, vRecoJet_tt2.Eta(), genJet->Phi, vRecoJet_tt2.Phi()) < MAX_MATCH_DIST) {
        iGenJet_tt2=iJet;
        genJet_tt2 = (Jet*) branchGenJet->At(iGenJet_tt2);
      }
      else if (deltaR(genJet->Eta, vRecoJet_6j1.Eta(), genJet->Phi, vRecoJet_6j1.Phi()) < MAX_MATCH_DIST) {
        iGenJet_6j1=iJet;
        genJet_6j1 = (Jet*) branchGenJet->At(iGenJet_6j1);
      }
      else if (deltaR(genJet->Eta, vRecoJet_6j2.Eta(), genJet->Phi, vRecoJet_6j2.Phi()) < MAX_MATCH_DIST) {
        iGenJet_6j2=iJet;
        genJet_6j2 = (Jet*) branchGenJet->At(iGenJet_6j2);
      }
    }

    if (genJetTau1) {
      ptTau1_genJet=genJetTau1->PT;
      etaTau1_genJet=genJetTau1->Eta;
      phiTau1_genJet=genJetTau1->Phi;
      mTau1_genJet=genJetTau1->Mass;
    }

    if (genJetTau2) {
      ptTau2_genJet=genJetTau2->PT;
      etaTau2_genJet=genJetTau2->Eta;
      phiTau2_genJet=genJetTau2->Phi;
      mTau2_genJet=genJetTau2->Mass;
    }

    if (genJet_tt1) {
      ptJet_tt1_gen=genJet_tt1->PT;
      etaJet_tt1_gen=genJet_tt1->Eta;
      phiJet_tt1_gen=genJet_tt1->Phi;
      mJet_tt1_gen=genJet_tt1->Mass;
    }

    if (genJet_tt2) {
      ptJet_tt2_gen=genJet_tt2->PT;
      etaJet_tt2_gen=genJet_tt2->Eta;
      phiJet_tt2_gen=genJet_tt2->Phi;
      mJet_tt2_gen=genJet_tt2->Mass;
    }

    if (genJet_6j1) {
      ptJet_6j1_gen=genJet_6j1->PT;
      etaJet_6j1_gen=genJet_6j1->Eta;
      phiJet_6j1_gen=genJet_6j1->Phi;
      mJet_6j1_gen=genJet_6j1->Mass;
    }

    if (genJet_6j2) {
      ptJet_6j2_gen=genJet_6j2->PT;
      etaJet_6j2_gen=genJet_6j2->Eta;
      phiJet_6j2_gen=genJet_6j2->Phi;
      mJet_6j2_gen=genJet_6j2->Mass;
    }


    if (0) { 
    //if ( (isBBBB==1 || isVBF4B==1) && ptH1_gen>0 && ptH2_gen>0 ) {

      LorentzVector vGenB1(ptB1_gen, etaB1_gen, phiB1_gen, mB1_gen);
      LorentzVector vGenB2(ptB2_gen, etaB2_gen, phiB2_gen, mB2_gen);
      LorentzVector vGenB3(ptB3_gen, etaB3_gen, phiB3_gen, mB3_gen);
      LorentzVector vGenB4(ptB4_gen, etaB4_gen, phiB4_gen, mB4_gen);
      
      LorentzVector vGenH1(ptH1_gen, etaH1_gen, phiH1_gen, mH1_gen);
      LorentzVector vGenH2(ptH2_gen, etaH2_gen, phiH2_gen, mH2_gen);
      
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
  cout << "Selection complete" << endl;

}

Float_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 ) {

  const Float_t pi = 3.14159265358979;

  Float_t etaDiff = (eta1-eta2);
  Float_t phiDiff = fabs(phi1-phi2);
  while ( phiDiff>pi ) phiDiff = fabs(phiDiff-2.0*pi);

  Float_t deltaRSquared = etaDiff*etaDiff + phiDiff*phiDiff;

  return TMath::Sqrt(deltaRSquared);

}

Float_t deltaPhi( const Float_t phi1, const Float_t phi2 ) {

  const Float_t pi = 3.14159265358979;

  Float_t phiDiff = fabs(phi1-phi2);
  while ( phiDiff>pi ) phiDiff = fabs(phiDiff-2.0*pi);

  return phiDiff;

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

double doJetcorr(mithep::jcorr *corrector,Jet* ijet,double rho_2,double rho_1,double rho_0)
{
  double area = TMath::Sqrt(ijet->AreaX*ijet->AreaX+ijet->AreaY*ijet->AreaY);
  //std::cout << "Area pT  "  << area  << std::endl;
  double rrho=0;
  if(fabs(ijet->Eta)>4.0) rrho=rho_2;
  else if(fabs(ijet->Eta)>2.5) rrho=rho_1;
  else rrho=rho_0;
  double corr = corrector->getCorrection(ijet->PT,ijet->Eta,rrho,area);
  //std::cout << ijet->PT << "  " <<  corr << std::endl;
  //return 1.0;
  return corr;
}
