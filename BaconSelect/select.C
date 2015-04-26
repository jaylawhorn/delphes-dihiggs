#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TClonesArray.h>           // ROOT array class
#include <TBenchmark.h>             // class to track macro running statistics
#include <TVector2.h>               // 2D vector class
#include <TMath.h>                  // ROOT math library
#include <TRandom3.h>
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <sstream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <TChain.h>
#include <TH1.h>
#include "Math/LorentzVector.h"     // 4-vector class

#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenJet.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"

#include "TauAnalysis/SVFitHelper/interface/TSVfit.h"
#include "TauAnalysis/SVFitHelper/interface/TSVfitter.h"
#include "mt2.hh"

#endif

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;

Float_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 );

Double_t jetRes(Double_t pt, Double_t eta);

Double_t muRes(Double_t pt, Double_t eta);

Double_t btag(Double_t pt, Double_t eta);

Double_t muEff(Double_t pt, Double_t eta);

void select(const TString //input="root://eoscms.cern.ch//store/group/upgrade/di_higgs_backgrounds/gFHHTobbtt_TuneZ2_14TeV_madgraph/Bacon/gfhhbbtautau_mad_1.root.root",
	    input="root://eoscms.cern.ch//store/group/upgrade/di_higgs_backgrounds/TT_jets/Bacon/tt_mad_9956.root",
	    const TString output="test.root") {

  TChain chain("Events");
  chain.Add(input);
  
  const Float_t DR=0.4;

  const Double_t MUON_MASS = 0.105658369;
  const Double_t ELE_MASS  = 0.000511;
  const Double_t TAU_MASS  = 1.77682;
  
  // tau decay modes
  enum { hadron=1, electron, muon };

  TRandom3 *rng = new TRandom3();
  rng->SetSeed(0);

  // setup mt2 minimizer
  ROOT::Math::Minimizer *min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
  min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
  min->SetTolerance(10.0);
  min->SetPrintLevel(0);

  const Float_t lcov00 = 2500;
  const Float_t lcov10 = 500;
  const Float_t lcov01 = 500;
  const Float_t lcov11 = 2500;
  const Float_t lcov00pp = 225;
  const Float_t lcov10pp = 45;
  const Float_t lcov01pp = 45;
  const Float_t lcov11pp = 225;

  TVector2 vtau1(0,0), vtau2(0,0), pumpt(0,0);
  TVector2 vb1(0,0), vb2(0,0);
  Double_t mt2pileup=0;

  // setup svfit
  mithep::TSVfitter *fitter = new mithep::TSVfitter();

  // Data structures to store info from TTrees
  baconhep::TGenEventInfo *info = new baconhep::TGenEventInfo();
  TClonesArray *jet             = new TClonesArray("baconhep::TGenJet");
  TClonesArray *part            = new TClonesArray("baconhep::TGenParticle");
  
  chain.SetBranchAddress("GenEvtInfo",  &info);        TBranch *infoBr     = chain.GetBranch("GenEvtInfo");
  chain.SetBranchAddress("GenJet" ,     &jet );        TBranch *jetBr      = chain.GetBranch("GenJet");
  chain.SetBranchAddress("GenParticle", &part);        TBranch *partBr     = chain.GetBranch("GenParticle");
  
  TFile *outfile = new TFile(output, "RECREATE");

  TTree *infoTree = new TTree("Count", "Count");
  Long64_t n = chain.GetEntries();
  infoTree->Branch("n", &n, "n/i");
  infoTree->Fill();

  TTree *outtree = new TTree("Events", "Events");

  Int_t eventType;
  Float_t eventWeight;
  Float_t metPx, metPy;
  Float_t met, metPhi;
  Float_t met2, met2Phi;
  Double_t m_sv;
  Double_t m_svpileup;

  Int_t tauCat1=0, tauCat2=0, tauIso1, tauIso2;
  Int_t nProngTau1=0, nProngTau2=0;
  Int_t tauDecay1=0, tauDecay2=0;

  Float_t ptTau1, ptTau2, ptB1, ptB2;
  Float_t etaTau1, etaTau2, etaB1, etaB2;
  Float_t phiTau1, phiTau2, phiB1, phiB2;
  Float_t mTau1, mTau2, mB1, mB2;

  Float_t ptTau1_gen, ptTau2_gen, ptB1_gen, ptB2_gen;
  Float_t etaTau1_gen, etaTau2_gen, etaB1_gen, etaB2_gen;
  Float_t phiTau1_gen, phiTau2_gen, phiB1_gen, phiB2_gen;
  Float_t mTau1_gen, mTau2_gen, mB1_gen, mB2_gen;

  Float_t mTT, mBB1, dRTT, dRBB1, mTT_gen, mBB1_gen;

  Int_t nBjetWithHadTau=0, nBjetWithMuon=0, nBjetWithEle=0;
  Int_t nLepHigherPt=0;

  Int_t nTauTau=0, nMuTau=0;

  //hadronic final states and "prong-ness"
  //decayModes.push_back(20213); isOneProng.push_back(-1); // a1(1260)
  //decayModes.push_back(213); isOneProng.push_back(1); // rho
  //decayModes.push_back(211); isOneProng.push_back(1); // pi+
  //decayModes.push_back(321); isOneProng.push_back(1); // K+
  //decayModes.push_back(323); isOneProng.push_back(1); // K*+
  
  outtree->Branch("eventWeight",    &eventWeight,    "eventWeight/f");  // event weight from cross-section and Event->Weight
  outtree->Branch("eventType",      &eventType,      "eventType/i");    // see line 47

  outtree->Branch("metPx",     &metPx,     "metPx/F");
  outtree->Branch("metPy",     &metPy,     "metPy/F");

  outtree->Branch("met",       &met,       "met/F");
  outtree->Branch("metPhi",    &metPhi,    "metPhi/F");
  outtree->Branch("met2",      &met2,      "met2/F");
  outtree->Branch("met2Phi",   &met2Phi,   "met2Phi/F");

  outtree->Branch("m_sv",           &m_sv,           "m_sv/D");         // "SVFit mass estimate" 
  outtree->Branch("m_svpileup",     &m_svpileup,     "m_svpileup/D");   // "SVFit mass estimate with pileup jet ID MET"

  outtree->Branch("ptTau1",         &ptTau1,         "ptTau1/f");       // pt(Tau1)
  outtree->Branch("etaTau1",        &etaTau1,        "etaTau1/f");      // eta(Tau1)
  outtree->Branch("phiTau1",        &phiTau1,        "phiTau1/f");      // phi(Tau1)
  outtree->Branch("mTau1",          &mTau1,          "mTau1/f");        // m(Tau1)

  outtree->Branch("ptTau2",         &ptTau2,         "ptTau2/f");       // pt(Tau2)
  outtree->Branch("etaTau2",        &etaTau2,        "etaTau2/f");      // eta(Tau2)
  outtree->Branch("phiTau2",        &phiTau2,        "phiTau2/f");      // phi(Tau2)
  outtree->Branch("mTau2",          &mTau2,          "mTau2/f");        // m(Tau2)

  outtree->Branch("ptB1",           &ptB1,           "ptB1/f");         // pt(B1)
  outtree->Branch("etaB1",          &etaB1,          "etaB1/f");        // eta(B1)
  outtree->Branch("phiB1",          &phiB1,          "phiB1/f");        // phi(B1)
  outtree->Branch("mB1",            &mB1,            "mB1/f");          // m(B1)

  outtree->Branch("ptB2",           &ptB2,           "ptB2/f");         // pt(B2)
  outtree->Branch("etaB2",          &etaB2,          "etaB2/f");        // eta(B2)
  outtree->Branch("phiB2",          &phiB2,          "phiB2/f");        // phi(B2)
  outtree->Branch("mB2",            &mB2,            "mB2/f");          // m(B2)

  outtree->Branch("ptTau1_gen",         &ptTau1_gen,         "ptTau1_gen/f");       // pt(Tau1)
  outtree->Branch("etaTau1_gen",        &etaTau1_gen,        "etaTau1_gen/f");      // eta(Tau1)
  outtree->Branch("phiTau1_gen",        &phiTau1_gen,        "phiTau1_gen/f");      // phi(Tau1)
  outtree->Branch("mTau1_gen",          &mTau1_gen,          "mTau1_gen/f");        // m(Tau1)
  outtree->Branch("tauCat1",            &tauCat1,            "tauCat1/i");      // leading tau final state - jet, muon, electron
  outtree->Branch("nProngTau1",         &nProngTau1,         "nProngTau1/i");
  outtree->Branch("tauDecay1",          &tauDecay1,          "tauDecay1/i");
  outtree->Branch("tauIso1",            &tauIso1,            "tauIso1/f");      // leading tau final state - jet, muon, electron

  outtree->Branch("ptTau2_gen",         &ptTau2_gen,         "ptTau2_gen/f");       // pt(Tau2)
  outtree->Branch("etaTau2_gen",        &etaTau2_gen,        "etaTau2_gen/f");      // eta(Tau2)
  outtree->Branch("phiTau2_gen",        &phiTau2_gen,        "phiTau2_gen/f");      // phi(Tau2)
  outtree->Branch("mTau2_gen",          &mTau2_gen,          "mTau2_gen/f");        // m(Tau2)
  outtree->Branch("tauCat2",            &tauCat2,            "tauCat2/i");      // second tau final state - jet, muon, electron
  outtree->Branch("nProngTau2",         &nProngTau2,         "nProngTau2/i");
  outtree->Branch("tauDecay2",          &tauDecay2,          "tauDecay2/i");
  outtree->Branch("tauIso2",            &tauIso2,            "tauIso2/f");      // second tau final state - jet, muon, electron

  outtree->Branch("ptB1_gen",           &ptB1_gen,           "ptB1_gen/f");         // pt(B1)
  outtree->Branch("etaB1_gen",          &etaB1_gen,          "etaB1_gen/f");        // eta(B1)
  outtree->Branch("phiB1_gen",          &phiB1_gen,          "phiB1_gen/f");        // phi(B1)
  outtree->Branch("mB1_gen",            &mB1_gen,            "mB1_gen/f");          // m(B1)

  outtree->Branch("ptB2_gen",           &ptB2_gen,           "ptB2_gen/f");         // pt(B2)
  outtree->Branch("etaB2_gen",          &etaB2_gen,          "etaB2_gen/f");        // eta(B2)
  outtree->Branch("phiB2_gen",          &phiB2_gen,          "phiB2_gen/f");        // phi(B2)
  outtree->Branch("mB2_gen",            &mB2_gen,            "mB2_gen/f");          // m(B2)

  outtree->Branch("mTT",                &mTT,                "mTT/f");          // m(TT)
  outtree->Branch("mTT_gen",            &mTT_gen,            "mTT_gen/f");          // m(TT)
  outtree->Branch("mBB1",               &mBB1,               "mBB1/f");         // m(BB1)
  outtree->Branch("mBB1_gen",           &mBB1_gen,           "mBB1_gen/f");         // m(BB1)
  outtree->Branch("mt2pileup",          &mt2pileup,          "mt2pileup/D");    // PUJetID "stransverse mass" (HH)

  outtree->Branch("dRBB1",              &dRBB1,              "dRBB1/f");
  outtree->Branch("dRTT",               &dRTT,               "dRTT/f");

  outtree->Branch("nBjetWithHadTau",    &nBjetWithHadTau,    "nBjetWithHadTau/i");
  outtree->Branch("nBjetWithMuon",      &nBjetWithMuon,      "nBjetWithMuon/i");
  outtree->Branch("nBjetWithEle",       &nBjetWithEle,       "nBjetWithEle/i");
  outtree->Branch("nLepHigherPt",       &nLepHigherPt,       "nLepHigherPt/i");

  Int_t nHadTau=0, nOneP=0, nThreeP=0;

  for (Int_t i=0; i<chain.GetEntries(); i++) {
  //for (Int_t i=15213; i<15213+1; i++) {
  //for (Int_t i=0; i<5000; i++) {
    infoBr->GetEntry(i);
    //cout << " ----" << endl;
    part->Clear(); partBr->GetEntry(i);
    jet->Clear(); jetBr->GetEntry(i);
    
    if (part->GetEntries()==0) continue;
    
    Int_t iTau1=-1, iTau2=-1, iEle1=-1, iEle2=-1, iMu1=-1, iMu2=-1;
    Int_t iB1=-1, iB2=-1, iT1=-1, iT2=-1;

    nBjetWithHadTau=0; nBjetWithMuon=0; nBjetWithEle=0;
    nLepHigherPt=0;

    metPx=-999; metPy=-999;
    met=-999; met2=-999; metPhi=-999; met2Phi=-999;
    tauCat1=0; tauCat2=0;
    nProngTau1=0; nProngTau2=0;

    ptTau1=-999; etaTau1=-999; phiTau1=-999; mTau1=-999;
    ptTau2=-999; etaTau2=-999; phiTau2=-999; mTau2=-999;
    ptB1=-999; etaB1=-999; phiB1=-999; mB1=-999;
    ptB2=-999; etaB2=-999; phiB2=-999; mB2=-999;

    ptTau1_gen=-999; etaTau1_gen=-999; phiTau1_gen=-999; mTau1_gen=-999; tauIso1=0;
    ptTau2_gen=-999; etaTau2_gen=-999; phiTau2_gen=-999; mTau2_gen=-999; tauIso2=0;
    ptB1_gen=-999; etaB1_gen=-999; phiB1_gen=-999; mB1_gen=-999;
    ptB2_gen=-999; etaB2_gen=-999; phiB2_gen=-999; mB2_gen=-999;
    mTT=-999; mBB1=-999; mt2pileup=-999; mTT_gen=-999; mBB1_gen=-999;
    eventType=2;

    eventWeight=1.24*341*1000; //fb^-1
    //eventWeight=2.92;

    metPx=info->metPx;
    metPy=info->metPy;

    TVector2 genMet(metPx, metPy);

    Float_t sMetPx=rng->Gaus(info->metPx, 20);
    Float_t sMetPy=rng->Gaus(info->metPy, 20);
    //cout << "----------" << endl;
    //cout << "metPx " << metPx << " " << sMetPx << endl;
    //cout << "metPy " << metPy << " " << sMetPy << endl;

    TVector2 smearedMet(sMetPx, sMetPy);

    met2=smearedMet.Mod();
    met2Phi=smearedMet.Phi();

    //cout << "Gmet " << genMet.Mod() << " " << genMet.Phi() << endl;
    //cout << "Smet " << smearedMet.Mod() << " " << smearedMet.Phi() << endl;

    Int_t isLep=0;
    for (Int_t j=0; j<part->GetEntries(); j++) {
      const baconhep::TGenParticle* genloop = (baconhep::TGenParticle*) ((*part)[j]);
      Int_t parentPdg=dynamic_cast<baconhep::TGenParticle *>(part->At(genloop->parent>-1 ? genloop->parent : 0))->pdgId;

      if (!(fabs(genloop->pdgId)==15)&&!(fabs(parentPdg)==15)) continue;
      if ( (fabs(genloop->pdgId)==13||fabs(genloop->pdgId)==11||fabs(genloop->pdgId)==24) && (fabs(parentPdg)==15) ) {
	isLep=1;
	continue;
      }

      if (fabs(genloop->pdgId)==15) {
	if (iTau1==-1) {
	  iTau1=j;
	}
	else if (genloop->parent==iTau1) {
	  iTau1=j;
	}
	else if (iTau2==-1) {
	  iTau2=j;
	}
	else if (genloop->parent==iTau2) {
	  iTau2=j;
	}
      }	
    }

    for (Int_t j=0; j<part->GetEntries(); j++) {
      const baconhep::TGenParticle* genloop = (baconhep::TGenParticle*) ((*part)[j]);
      Int_t parentPdg=dynamic_cast<baconhep::TGenParticle *>(part->At(genloop->parent>-1 ? genloop->parent : 0))->pdgId;
      
      if (fabs(genloop->pdgId)==13 && genloop->status==1) {
	if (iMu1==-1) iMu1=j;
	else if (genloop->pt>dynamic_cast<baconhep::TGenParticle *>(part->At(iMu1))->pt) iMu1=j;
      }
    }
    
    /*    if (iTau1>0 && iTau2>0 && iMu1>0) {
      //removing tau particles that decayed to muons
      Int_t iTemp1=iTau1, iTemp2=iTau2;
      const baconhep::TGenParticle* mu = (baconhep::TGenParticle*) ((*part)[iTemp1]); 
      while (iTemp1>0) {
	//cout << iTemp1 << endl;
	const baconhep::TGenParticle* genloop = (baconhep::TGenParticle*) ((*part)[iTemp1]); 
	if (iTemp1 == mu->parent) { iTemp1=0; iTau1=-1; }
	else { iTemp1=genloop->parent; }
      }
      //cout << "." << endl;
      while (iTemp2>0) {
	//cout << iTemp2 << endl;
	const baconhep::TGenParticle* genloop = (baconhep::TGenParticle*) ((*part)[iTemp2]); 
	if (iTemp2 == mu->parent) { iTemp2=0; iTau2=-1; }
	else { iTemp2=genloop->parent; }
      }      
      }*/
    
    //if (iTau1==-1 && iTau2>0) {iTau1=iTau2; iTau2=-1; }

    if (iTau1==-1 || (iMu1==-1&&iTau2==-1)) continue;
    if (isLep>0 && iMu1==-1) continue;
    //else cout << "wtf? " << iTau1 << " " << iTau2 << " " << iMu1 << endl;

    /*    if (iTau1>0 && iTau2>0 && isLep==0) {
      cout << i << " found hadronic taus at " << iTau1 << " and " << iTau2 << endl;
    }
    if (iTau1>0 && iMu1>0) {
      cout << i << " found a hadronic tau at " << iTau1 << " and muon at " << iMu1 << endl;
      }*/

    Int_t aPos1=-1, aPos2=-1;

    for (Int_t j=0; j<part->GetEntries(); j++) {
      const baconhep::TGenParticle* genloop = (baconhep::TGenParticle*) ((*part)[j]);
      
      Int_t pdg=fabs(genloop->pdgId);
      Int_t parentPdg=fabs(dynamic_cast<baconhep::TGenParticle *>(part->At(genloop->parent>-1 ? genloop->parent : 0))->pdgId);

      if (genloop->parent==iTau1 && pdg>22) {
	//cout << "found a charged hadronic for tau 1: " << pdg << endl;
	tauDecay1=pdg;
	if (pdg==20213) aPos1=j;
	else nProngTau1++;
      }
      else if ((aPos1>0) && (genloop->parent==aPos1) && pdg!=111 && pdg>22) {
	//cout << "found an a1 decay: " << pdg << endl;
	nProngTau1++;
      }
      else if (iTau2>-1 && genloop->parent==iTau2 && pdg>22) {
	tauDecay2=pdg;
	//cout << "found a charged hadronic for tau 2: " << pdg << endl;
	if (pdg==20213) aPos2=j;
	else nProngTau2++;
      }
      else if (iTau2>-1 && (aPos2>0) && (genloop->parent==aPos2) && pdg!=111 && pdg>22) {
	//cout << "found an a1 decay: " << pdg << endl;
	nProngTau2++;
      }      
    }

    //cout << nProngTau1 << " " << nProngTau2 << endl;

    if (nProngTau1==1) nOneP++;
    else if (nProngTau1==3) nThreeP++;

    if (nProngTau2==1) nOneP++;
    else if (nProngTau2==3) nThreeP++;
    
    const baconhep::TGenParticle* genTau1 = (baconhep::TGenParticle*) ((*part)[iTau1]);
    LorentzVector vGenTau1(genTau1->pt,genTau1->eta,genTau1->phi,genTau1->mass);
    
    const baconhep::TGenParticle* genTau2 = (baconhep::TGenParticle*) ((*part)[ (iTau2>0 ? iTau2 : iMu1) ]);
    LorentzVector vGenTau2(genTau2->pt,genTau2->eta,genTau2->phi,genTau2->mass);
    
    //cout << vGenTau1.Pt() << " " << vGenTau1.Eta() << " " << vGenTau1.Phi() << " " << vGenTau1.M() << endl;
    //cout << vGenTau2.Pt() << " " << vGenTau2.Eta() << " " << vGenTau2.Phi() << " " << vGenTau2.M() << " " << iMu1 << endl;
    
    if (deltaR(genTau1->eta, genTau2->eta, genTau1->phi, genTau2->phi)<0.4) continue;

    for (Int_t j=0; j<part->GetEntries(); j++) {
      const baconhep::TGenParticle* genloop = (baconhep::TGenParticle*) ((*part)[j]);
      if (fabs(genloop->pdgId)!=16) continue;
      LorentzVector vTemp(genloop->pt, genloop->eta, genloop->phi, 0);
      if (genloop->parent==iTau1) {
        vGenTau1=vGenTau1-vTemp;
      }
      if (genloop->parent==iTau2) {
	vGenTau2=vGenTau2-vTemp;
      }
    }

    //cout << vGenTau1.Pt() << " " << vGenTau1.Eta() << " " << vGenTau1.Phi() << " " << vGenTau1.M() << endl;
    //cout << vGenTau2.Pt() << " " << vGenTau2.Eta() << " " << vGenTau2.Phi() << " " << vGenTau2.M() << " " << iMu1 << endl;

    for (Int_t j=0; j<jet->GetEntries(); j++) {
      const baconhep::TGenJet* loop = (baconhep::TGenJet*) ((*jet)[j]);

      if (iTau1>-1) {
	if (deltaR(loop->eta,vGenTau1.Eta(),loop->phi,vGenTau1.Phi())<0.4) iT1=j;
      }
      if (iTau2>-1) {
	if (deltaR(loop->eta,vGenTau2.Eta(),loop->phi,vGenTau2.Phi())<0.4) iT2=j;
      }
    }

    //cout << iTau1 << " " << iTau2 << " " << iMu1 << ", " << iT1 << " " << iT2 << endl;

    if (iT1==-1||(iT2==-1&&iMu1==-1)) continue;
    if (isLep>0 && iMu1==-1) continue;
    //if (iT2==-1) continue;
    
    //if (iT1>-1 && iT2>-1 && isLep==0) {
    //nTauTau++;
      //cout << i << ": " << iTau1 << " " << iTau2 << ", " << iT1 << " " << iT2 << endl;
    //}
    //else if (iT1>-1 && iMu1>-1) {
    //nMuTau++;
    //}
    //else cout << iT1 << " " << iT2 << " " << iMu1 << endl;

    LorentzVector genSumPt(0,0,0,0);
    LorentzVector smeSumPt(0,0,0,0);

    for (Int_t j=0; j<jet->GetEntries(); j++) {
      const baconhep::TGenJet* loop = (baconhep::TGenJet*) ((*jet)[j]);
      LorentzVector vTempG(loop->pt, loop->eta, loop->phi, loop->mass);
      LorentzVector vTempS(loop->pt*(rng->Gaus(1,jetRes(loop->pt, loop->eta))), loop->eta, loop->phi, loop->mass);
      genSumPt+=vTempG;
      smeSumPt+=vTempS;
      
      if (fabs(loop->mcFlavor)==5 && iB1==-1) iB1=j;
      else if (fabs(loop->mcFlavor)==5 && iB2==-1) iB2=j;
    }

    if (iB1>-1) {
      const baconhep::TGenJet *bjet = (baconhep::TGenJet*) ((*jet)[iB1]);
      ptB1_gen=bjet->pt;
      ptB1=ptB1_gen*(rng->Gaus(1,jetRes(bjet->pt, bjet->eta)));
      etaB1_gen=bjet->eta;
      etaB1=bjet->eta;
      mB1_gen=bjet->mass;
      mB1=bjet->mass;
      phiB1_gen=bjet->phi;
      phiB1=bjet->phi;
      eventWeight*=btag(ptB1, etaB1);
    }
    if (iB2>-1) {
      const baconhep::TGenJet *bjet = (baconhep::TGenJet*) ((*jet)[iB2]);
      ptB2_gen=bjet->pt;
      ptB2=ptB2_gen*(rng->Gaus(1,jetRes(bjet->pt, bjet->eta)));
      etaB2_gen=bjet->eta;
      etaB2=bjet->eta;
      mB2_gen=bjet->mass;
      mB2=bjet->mass;
      phiB2_gen=bjet->phi;
      phiB2=bjet->phi;
      eventWeight*=btag(ptB2, etaB2);
    }

    if (iT1>-1 && iT2>-1 && isLep==0) {
      nTauTau++;
      const baconhep::TGenJet *taujet = (baconhep::TGenJet*) ((*jet)[iT1]);
      const baconhep::TGenJet *taujet2 = (baconhep::TGenJet*) ((*jet)[iT2]);
      
      Float_t fact1=(rng->Gaus(1,jetRes(taujet->pt, taujet->eta)));
      Float_t fact2=(rng->Gaus(1,jetRes(taujet2->pt, taujet2->eta)));
      
      Float_t ptTauJet1=taujet->pt*fact1;
      Float_t ptTauJet2=taujet2->pt*fact2;
      
      if (ptTauJet1>ptTauJet2) {
	ptTau1_gen=taujet->pt;
	ptTau1=ptTauJet1;
	etaTau1_gen=taujet->eta;
	etaTau1=taujet->eta;
	mTau1_gen=taujet->mass;
	mTau1=taujet->mass;
	phiTau1_gen=taujet->phi;
	phiTau1=taujet->phi;
	tauCat1=1;
	
	ptTau2_gen=taujet2->pt;
	ptTau2=ptTauJet2;
	etaTau2_gen=taujet2->eta;
	etaTau2=taujet2->eta;
	mTau2_gen=taujet2->mass;
	mTau2=taujet2->mass;
	phiTau2_gen=taujet2->phi;
	phiTau2=taujet2->phi;
	tauCat2=1;
      }
      else {
	ptTau2_gen=taujet2->pt;
	ptTau2=ptTauJet1;
	etaTau2_gen=taujet->eta;
	etaTau2=taujet->eta;
	mTau2_gen=taujet->mass;
	mTau2=taujet->mass;
	phiTau2_gen=taujet->phi;
	phiTau2=taujet->phi;
	tauCat2=1;
	
	ptTau1_gen=taujet->pt;
	ptTau1=ptTauJet2;
	etaTau1_gen=taujet2->eta;
	etaTau1=taujet2->eta;
	mTau1_gen=taujet2->mass;
	mTau1=taujet2->mass;
	phiTau1_gen=taujet2->phi;
	phiTau1=taujet2->phi;
	tauCat1=1;
      }

      if (deltaR(etaB1,taujet->eta,phiB1,taujet->phi)<0.4||deltaR(etaB1,taujet2->eta,phiB1,taujet2->phi)<0.4) nBjetWithHadTau++;
      if (deltaR(etaB2,taujet->eta,phiB2,taujet->phi)<0.4||deltaR(etaB2,taujet2->eta,phiB2,taujet2->phi)<0.4) nBjetWithHadTau++;

      for (Int_t j=0; j<part->GetEntries(); j++) {
	const baconhep::TGenParticle* genloop = (baconhep::TGenParticle*) ((*part)[j]);
	if ( genloop->pt<20) continue;
	if ( fabs(genloop->pdgId)!=13&&fabs(genloop->pdgId)!=11 ) continue;
	
	if ((fabs(genloop->pdgId)==11) && (deltaR(etaB1,genloop->eta,phiB1,genloop->phi)<0.4||deltaR(etaB1,genloop->eta,phiB1,genloop->phi)<0.4)) nBjetWithEle++;
	if ((fabs(genloop->pdgId)==13) && (deltaR(etaB1,genloop->eta,phiB1,genloop->phi)<0.4||deltaR(etaB1,genloop->eta,phiB1,genloop->phi)<0.4)) nBjetWithMuon++;
	
	if (genloop->pt>TMath::Min(ptTau1,ptTau2))nLepHigherPt++;
      }
      
      eventWeight*=0.65*0.65;
    }
    else if (iT1>-1 && iMu1>-1) {
      nMuTau++;
      const baconhep::TGenJet *taujet = (baconhep::TGenJet*) ((*jet)[iT1]);
      const baconhep::TGenParticle* mutau = (baconhep::TGenParticle*) ((*part)[iMu1]);

      Float_t fact1=(rng->Gaus(1,jetRes(taujet->pt, taujet->eta)));
      Float_t fact2=(rng->Gaus(1,muRes(mutau->pt, mutau->eta)));

      Float_t ptTauJet1=taujet->pt*fact1;
      Float_t ptMuTau=mutau->pt*fact2;
      
      ptTau1_gen=taujet->pt;
      ptTau1=ptTauJet1;
      etaTau1_gen=taujet->eta;
      etaTau1=taujet->eta;
      mTau1_gen=taujet->mass;
      mTau1=taujet->mass;
      phiTau1_gen=taujet->phi;
      phiTau1=taujet->phi;
      tauCat1=1;

      ptTau2_gen=mutau->pt;
      ptTau2=ptMuTau;
      etaTau2_gen=mutau->eta;
      etaTau2=mutau->eta;
      mTau2_gen=mutau->mass;
      mTau2=mutau->mass;
      phiTau2_gen=mutau->phi;
      phiTau2=mutau->phi;
      tauCat2=3;
      
      if (deltaR(etaB1,taujet->eta,phiB1,taujet->phi)<0.4||deltaR(etaB2,taujet->eta,phiB2,taujet->phi)<0.4) nBjetWithHadTau++;

      for (Int_t j=0; j<part->GetEntries(); j++) {
	const baconhep::TGenParticle* genloop = (baconhep::TGenParticle*) ((*part)[j]);
	if ( genloop->pt<20) continue;
	if ( fabs(genloop->pdgId)!=13&&fabs(genloop->pdgId)!=11 ) continue;
	
	if ((fabs(genloop->pdgId)==11) && (deltaR(etaB1,genloop->eta,phiB1,genloop->phi)<0.4||deltaR(etaB1,genloop->eta,phiB1,genloop->phi)<0.4)) nBjetWithEle++;
	if ((fabs(genloop->pdgId)==13) && (deltaR(etaB1,genloop->eta,phiB1,genloop->phi)<0.4||deltaR(etaB1,genloop->eta,phiB1,genloop->phi)<0.4)) nBjetWithMuon++;
	
	if (genloop->pt>TMath::Min(ptTau1,ptTau2))nLepHigherPt++;
      }
      
      eventWeight*=muEff(ptTau2, etaTau2)*0.65;

    }

    if (ptB1==-999 || ptB2==-999 || ptTau1==-999 || ptTau2==-999) continue;

    //cout << " smeared sum " << endl;
    //cout << smeSumPt.Pt() << " " << smeSumPt.Eta() << " " << smeSumPt.Phi() << endl;
    //cout << " unsmeared sum " << endl;
    //cout << genSumPt.Pt() << " " << genSumPt.Eta() << " " << genSumPt.Phi() << endl;

    LorentzVector metCorr = smeSumPt - genSumPt;
    TVector2 metCorrT(metCorr.Pt(), metCorr.Phi());
    TVector2 newMet=genMet-metCorrT;
    met=newMet.Mod();
    metPhi=newMet.Phi();

    //cout << "MET? " << met << " " << metPhi << endl;

    LorentzVector recoB1(ptB1, etaB1, phiB1, mB1);
    LorentzVector recoB2(ptB2, etaB2, phiB2, mB2);
    LorentzVector recoTau1(ptTau1, etaTau1, phiTau1, mTau1);
    LorentzVector recoTau2(ptTau2, etaTau2, phiTau2, mTau2);

    LorentzVector vgenB1(ptB1_gen, etaB1, phiB1, mB1);
    LorentzVector vgenB2(ptB2_gen, etaB2, phiB2, mB2);
    LorentzVector vgenTau1(ptTau1_gen, etaTau1, phiTau1, mTau1);
    LorentzVector vgenTau2(ptTau2_gen, etaTau2, phiTau2, mTau2);

    dRTT=deltaR(recoTau1.Eta(), recoTau2.Eta(), recoTau1.Phi(), recoTau2.Phi());
    dRBB1=deltaR(recoB1.Eta(), recoB2.Eta(), recoB1.Phi(), recoB2.Phi());

    LorentzVector vTT_gen=vgenTau1+vgenTau2;
    mTT_gen=vTT_gen.M();

    LorentzVector vTT=recoTau1+recoTau2;
    mTT=vTT.M();

    LorentzVector vBB_gen=vgenB1+vgenB2;
    mBB1_gen=vBB_gen.M();

    LorentzVector vBB=recoB1+recoB2;
    mBB1=vBB.M();

    //if (0) {
    if (recoTau1.Pt()>0 && recoTau2.Pt()>0 && recoB1.Pt()>0 && recoB2.Pt()>0) {
      vtau1.SetMagPhi(recoTau1.Pt(), recoTau1.Phi());
      vtau2.SetMagPhi(recoTau2.Pt(), recoTau2.Phi());
      
      vb1.SetMagPhi(recoB1.Pt(), recoB1.Phi());
      vb2.SetMagPhi(recoB2.Pt(), recoB2.Phi());
      
      pumpt.Set(metPx, metPy);
      TVector2 puSumPt = vtau1+vtau2+pumpt;
      
      smT2 calcmt2 = smT2();
      calcmt2.SetB1(vb1);
      calcmt2.SetB2(vb2);
      calcmt2.SetMPT(puSumPt);
      calcmt2.SetMB1(mB1);
      calcmt2.SetMB2(mB2);
      calcmt2.SetMT1(mTau1);
      calcmt2.SetMT2(mTau2);
      
      TVector2 c1=puSumPt;
      TVector2 c2=puSumPt-c1;
      
      double step[2] = {0.1, 0.1};
      double variable[2] = { 0.5*c1.Mod(), 0.0 };
      
      ROOT::Math::Functor f(calcmt2,2);
      
      min->SetFunction(f);
      min->SetLimitedVariable(0,"cT",variable[0], step[0], 0.0, puSumPt.Mod());
      min->SetLimitedVariable(1,"cPhi",variable[1], step[1], 0.0, TMath::Pi());
      
      min->Minimize();
      mt2pileup = min->MinValue();
    }

    // ***********************************
    // Let's start with SVFit calculations
    // ***********************************
    
    //if (0) {
    if (recoTau1.Pt()>0 && recoTau2.Pt()>0) {
      int channel=0;
      mithep::TSVfit svfit;
      svfit.cov_00=lcov00pp;
      svfit.cov_01=lcov01pp;
      svfit.cov_10=lcov10pp;
      svfit.cov_11=lcov11pp;
      TLorentzVector lvec1;
      lvec1.SetPtEtaPhiM(ptTau1,etaTau1,phiTau1,mTau1);
      TLorentzVector lvec2; 
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
      //cout << "met " << met << " " << metPhi << endl;
      //m_sv = fitter->integrate(&svfit,met,metPhi,channel);
      //svfit.cov_00=lcov00pp;
      //svfit.cov_01=lcov01pp;
      //svfit.cov_10=lcov10pp;
      //svfit.cov_11=lcov11pp;
      //cout << "met2 " << met2 << " " << met2Phi << endl;
      m_svpileup = fitter->integrate(&svfit,met2,met2Phi,channel);
      //std::cout << tauCat1 << "  " <<  tauCat2 << "   " << channel  << "  "  << m_sv << "  "  <<  m_svpileup << "  " << std::endl; 
    }

    outtree->Fill();
    
  }
  
  //cout << chain.GetEntries() << " " << nTauTau << " " << nMuTau << endl;
  
  //cout << "total had. taus: " << nHadTau << endl;
  //cout << "one prong:       " << nOneP << " (" << nOneP*100/nHadTau << "%)" << endl;
  //cout << "three prong:     " << nThreeP << " (" << nThreeP*100/nHadTau << "%)" << endl;
  
  outfile->Write();
  outfile->Save();
  
}

Float_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 ) {

  const Float_t pi = 3.14159265358979;

  Float_t etaDiff = (eta1-eta2);
  Float_t phiDiff = fabs(phi1-phi2);
  while ( phiDiff>pi ) phiDiff = fabs(phiDiff-2.0*pi);

  Float_t deltaRSquared = etaDiff*etaDiff + phiDiff*phiDiff;

  return TMath::Sqrt(deltaRSquared);

}

Double_t jetRes(Double_t pt, Double_t eta) {

  if (fabs(eta)<3.0) return 0.15;
  else if (fabs(eta)<5.0) return 0.3;
  else return 0;

}

Double_t muRes(Double_t pt, Double_t eta) {

  return ( (abs(eta) <= 1.5) * (pt > 0.1   && pt <= 1.0)   * (0.015) + 
	   (abs(eta) <= 1.5) * (pt > 1.0   && pt <= 1.0e1) * (0.012) + 
	   (abs(eta) <= 1.5) * (pt > 1.0e1 && pt <= 2.0e2) * (0.015) + 
	   (abs(eta) <= 1.5) * (pt > 2.0e2)                * (0.03) + 
	   (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 0.1   && pt <= 1.0)   * (0.015) + 
	   (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0   && pt <= 1.0e1) * (0.015) + 
	   (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0e1 && pt <= 2.0e2) * (0.025) + 
	   (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 2.0e2)                * (0.03) +
	   (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 0.1   && pt <= 1.0)   * (0.017) + 
	   (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 1.0   && pt <= 10.0)  * (0.03) + 
	   (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 10.0  && pt <= 100.0) * (0.05) + 
	   (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 100.0)                * (0.30) + 
	   (abs(eta) > 3.0 && abs(eta) <= 3.5) * (pt > 0.1   && pt <= 1.0)   * (0.02) + 
	   (abs(eta) > 3.0 && abs(eta) <= 3.5) * (pt > 1.0   && pt <= 10.0)  * (0.04) + 
	   (abs(eta) > 3.0 && abs(eta) <= 3.5) * (pt > 10.0  && pt <= 100.0) * (0.07) + 
	   (abs(eta) > 3.0 && abs(eta) <= 3.5) * (pt > 100.0)                * (0.30) + 
	   (abs(eta) > 3.5 && abs(eta) <= 4.0) * (pt > 0.1   && pt <= 1.0)   * (0.025) + 
	   (abs(eta) > 3.5 && abs(eta) <= 4.0) * (pt > 1.0   && pt <= 10.0)  * (0.05) + 
	   (abs(eta) > 3.5 && abs(eta) <= 4.0) * (pt > 10.0  && pt <= 100.0) * (0.20) + 
	   (abs(eta) > 3.5 && abs(eta) <= 4.0) * (pt > 100.0)                * (0.80) );
}

Double_t btag(Double_t pt, Double_t eta) {

  return ( (pt <= 20.0) * (0.000) +
           (abs(eta) <= 1.8) * (pt > 20.0 && pt <= 30) * (0.536) +
           (abs(eta) <= 1.8) * (pt > 30.0 && pt <= 40) * (0.6439) +
           (abs(eta) <= 1.8) * (pt > 40.0 && pt <= 50) * (0.6504) +
           (abs(eta) <= 1.8) * (pt > 50.0 && pt <= 60) * (0.6716) +
           (abs(eta) <= 1.8) * (pt > 60.0 && pt <= 70) * (0.6841) +
           (abs(eta) <= 1.8) * (pt > 70.0 && pt <= 80) * (0.6896) +
           (abs(eta) <= 1.8) * (pt > 80.0 && pt <= 90) * (0.6916) +
           (abs(eta) <= 1.8) * (pt > 90.0 && pt <= 100) * (0.6882) +
           (abs(eta) <= 1.8) * (pt > 100.0 && pt <= 120) * (0.6838) +
           (abs(eta) <= 1.8) * (pt > 120.0 && pt <= 140) * (0.6715) +
           (abs(eta) <= 1.8) * (pt > 140.0 && pt <= 160) * (0.6554) +
           (abs(eta) <= 1.8) * (pt > 160.0 && pt <= 180) * (0.6366) +
           (abs(eta) <= 1.8) * (pt > 180.0 && pt <= 200) * (0.6192) +
           (abs(eta) <= 1.8) * (pt > 200.0 && pt <= 250) * (0.595) +
           (abs(eta) <= 1.8) * (pt > 250.0 && pt <= 300) * (0.5551) +
           (abs(eta) <= 1.8) * (pt > 300.0 && pt <= 350) * (0.5138) +
           (abs(eta) <= 1.8) * (pt > 350.0 && pt <= 400) * (0.4884) +
           (abs(eta) <= 1.8) * (pt > 400.0 && pt <= 500) * (0.4009) +
           (abs(eta) <= 1.8) * (pt > 500.0 && pt <= 600) * (0.3459) +
           (abs(eta) <= 1.8) * (pt > 600.0 && pt <= 700) * (0.2523) +
           (abs(eta) <= 1.8) * (pt > 700.0 && pt <= 800) * (0.2404) +
           (abs(eta) <= 1.8) * (pt > 800.0 && pt <= 1000) * (0.2198) +
           (abs(eta) <= 1.8) * (pt > 1000.0 && pt <= 1400) * (0.2263) +
           (abs(eta) <= 1.8) * (pt > 1400.0 && pt <= 2000) * (0.2614) +
           (abs(eta) <= 1.8) * (pt > 2000.0) * (0.3194) +
           (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt <= 20.0) * (0.000) +
           (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 20.0 && pt <= 30) * (0.3254) +
           (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 30.0 && pt <= 40) * (0.4339) +
           (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 40.0 && pt <= 50) * (0.4499) +
           (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 50.0 && pt <= 60) * (0.4716) +
           (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 60.0 && pt <= 70) * (0.4766) +
           (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 70.0 && pt <= 80) * (0.4788) +
           (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 80.0 && pt <= 90) * (0.4863) +
           (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 90.0 && pt <= 100) * (0.4891) +
           (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 100.0 && pt <= 120) * (0.462) +
           (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 120.0 && pt <= 140) * (0.4583) +
           (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 140.0 && pt <= 160) * (0.4247) +
           (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 160.0 && pt <= 180) * (0.3775) +
           (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 180.0 && pt <= 200) * (0.3734) +
           (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 200.0 && pt <= 250) * (0.3348) +
           (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 250.0 && pt <= 300) * (0.2939) +
           (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 300.0 && pt <= 350) * (0.285) +
           (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 350.0 && pt <= 400) * (0.2421) +
           (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 400.0 && pt <= 500) * (0.1565) +
           (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 500.0 && pt <= 600) * (0.1522) +
           (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 600.0 && pt <= 700) * (0.1231) +
           (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 700.0 && pt <= 800) * (0.1607) +
           (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 800.0 && pt <= 1000) * (0.1323) +
           (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1000.0 && pt <= 1400) * (0.0) +
           (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1400.0 && pt <= 2000) * (0.0) +
           (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 2000.0) * (0.0) +
           (abs(eta) > 2.4) * (0.000));
}

Double_t muEff(Double_t pt, Double_t eta) {

  return ((pt <= 0.2) * (0.00) + 
	  (abs(eta) <= 1.2) * (pt > 0.2 && pt <= 1.0) * (pt * 0.998) + 
	  (abs(eta) <= 1.2) * (pt > 1.0) * (0.998) + 
	  (abs(eta) > 1.2 && abs(eta) <= 2.5) * (pt > 0.2 && pt <= 1.0) * (pt*0.99) + 
	  (abs(eta) > 1.2 && abs(eta) <= 2.5) * (pt > 1.0) * (0.99) + 
	  (abs(eta) > 2.5 && abs(eta) <= 4.0) * (pt > 0.2 && pt <= 1.0) * (pt*0.95) + 
	  (abs(eta) > 2.5 && abs(eta) <= 4.0) * (pt > 1.0) * (0.95) + 
	  (abs(eta) > 4.0) * (0.00));

}
