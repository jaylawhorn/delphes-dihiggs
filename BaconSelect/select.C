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

Double_t res(Double_t pt, Double_t eta);

Double_t btag(Double_t pt, Double_t eta);

void select(const TString input="root://eoscms.cern.ch//store/group/upgrade/di_higgs_backgrounds/gFHHTobbtt_TuneZ2_14TeV_madgraph/Bacon/gfhhbbtautau_mad_90001.root.root", 
	    //input="root://eoscms.cern.ch//store/group/upgrade/di_higgs_backgrounds/TT_jets/Bacon/tt_mad_950.root",
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
  
  //cout << chain.GetEntries() << endl;

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

  Float_t ptTau1, ptTau2, ptB1, ptB2, ptC1, ptC2, ptJ1, ptJ2;
  Float_t etaTau1, etaTau2, etaB1, etaB2, etaC1, etaC2, etaJ1, etaJ2;
  Float_t phiTau1, phiTau2, phiB1, phiB2, phiC1, phiC2, phiJ1, phiJ2;
  Float_t mTau1, mTau2, mB1, mB2, mC1, mC2, mJ1, mJ2;

  Float_t ptTau1_gen, ptTau2_gen, ptB1_gen, ptB2_gen, ptC1_gen, ptC2_gen, ptJ1_gen, ptJ2_gen;
  Float_t etaTau1_gen, etaTau2_gen, etaB1_gen, etaB2_gen, etaC1_gen, etaC2_gen, etaJ1_gen, etaJ2_gen;
  Float_t phiTau1_gen, phiTau2_gen, phiB1_gen, phiB2_gen, phiC1_gen, phiC2_gen, phiJ1_gen, phiJ2_gen;
  Float_t mTau1_gen, mTau2_gen, mB1_gen, mB2_gen, mC1_gen, mC2_gen, mJ1_gen, mJ2_gen;

  Float_t mTT, mBB1;

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
  outtree->Branch("tauIso1",            &tauIso1,            "tauIso1/f");      // leading tau final state - jet, muon, electron

  outtree->Branch("ptTau2_gen",         &ptTau2_gen,         "ptTau2_gen/f");       // pt(Tau2)
  outtree->Branch("etaTau2_gen",        &etaTau2_gen,        "etaTau2_gen/f");      // eta(Tau2)
  outtree->Branch("phiTau2_gen",        &phiTau2_gen,        "phiTau2_gen/f");      // phi(Tau2)
  outtree->Branch("mTau2_gen",          &mTau2_gen,          "mTau2_gen/f");        // m(Tau2)
  outtree->Branch("tauCat2",            &tauCat2,            "tauCat2/i");      // second tau final state - jet, muon, electron
  outtree->Branch("tauIso2",            &tauIso2,            "tauIso2/f");      // second tau final state - jet, muon, electron

  outtree->Branch("ptB1_gen",           &ptB1_gen,           "ptB1_gen/f");         // pt(B1)
  outtree->Branch("etaB1_gen",          &etaB1_gen,          "etaB1_gen/f");        // eta(B1)
  outtree->Branch("phiB1_gen",          &phiB1_gen,          "phiB1_gen/f");        // phi(B1)
  outtree->Branch("mB1_gen",            &mB1_gen,            "mB1_gen/f");          // m(B1)

  outtree->Branch("ptB2_gen",           &ptB2_gen,           "ptB2_gen/f");         // pt(B2)
  outtree->Branch("etaB2_gen",          &etaB2_gen,          "etaB2_gen/f");        // eta(B2)
  outtree->Branch("phiB2_gen",          &phiB2_gen,          "phiB2_gen/f");        // phi(B2)
  outtree->Branch("mB2_gen",            &mB2_gen,            "mB2_gen/f");          // m(B2)

  outtree->Branch("ptC1_gen",           &ptC1_gen,           "ptC1_gen/f");         // pt(C1)
  outtree->Branch("etaC1_gen",          &etaC1_gen,          "etaC1_gen/f");        // eta(C1)
  outtree->Branch("phiC1_gen",          &phiC1_gen,          "phiC1_gen/f");        // phi(C1)
  outtree->Branch("mC1_gen",            &mC1_gen,            "mC1_gen/f");          // m(C1)

  outtree->Branch("ptC2_gen",           &ptC2_gen,           "ptC2_gen/f");         // pt(C2)
  outtree->Branch("etaC2_gen",          &etaC2_gen,          "etaC2_gen/f");        // eta(C2)
  outtree->Branch("phiC2_gen",          &phiC2_gen,          "phiC2_gen/f");        // phi(C2)
  outtree->Branch("mC2_gen",            &mC2_gen,            "mC2_gen/f");          // m(C2)

  outtree->Branch("ptJ1_gen",           &ptJ1_gen,           "ptJ1_gen/f");         // pt(J1)
  outtree->Branch("etaJ1_gen",          &etaJ1_gen,          "etaJ1_gen/f");        // eta(J1)
  outtree->Branch("phiJ1_gen",          &phiJ1_gen,          "phiJ1_gen/f");        // phi(J1)
  outtree->Branch("mJ1_gen",            &mJ1_gen,            "mJ1_gen/f");          // m(J1)

  outtree->Branch("ptJ2_gen",           &ptJ2_gen,           "ptJ2_gen/f");         // pt(J2)
  outtree->Branch("etaJ2_gen",          &etaJ2_gen,          "etaJ2_gen/f");        // eta(J2)
  outtree->Branch("phiJ2_gen",          &phiJ2_gen,          "phiJ2_gen/f");        // phi(J2)
  outtree->Branch("mJ2_gen",            &mJ2_gen,            "mJ2_gen/f");          // m(J2)

  outtree->Branch("mTT",            &mTT,            "mTT/f");          // m(TT)
  outtree->Branch("mBB1",           &mBB1,           "mBB1/f");         // m(BB1)
  outtree->Branch("mt2pileup",      &mt2pileup,      "mt2pileup/D");    // PUJetID "stransverse mass" (HH)

  Int_t iHH=0, iEH=0, iMH=0;
  
  for (Int_t i=0; i<chain.GetEntries(); i++) {
  //for (Int_t i=0; i<100; i++) {
    infoBr->GetEntry(i);
    
    part->Clear(); partBr->GetEntry(i);
    jet->Clear(); jetBr->GetEntry(i);
    
    if (part->GetEntries()==0) continue;
    
    Int_t iTau1=-1, iTau2=-1, iEle1=-1, iEle2=-1, iMu1=-1, iMu2=-1;
    Int_t iB1=-1, iB2=-1, iC1=-1, iC2=-1, iT1=-1, iT2=-1, iJ1=-1, iJ2=-1;

    metPx=-999; metPy=-999;
    met=-999; met2=-999; metPhi=-999; met2Phi=-999;
    tauCat1=0; tauCat2=0;

    ptTau1=-999; etaTau1=-999; phiTau1=-999; mTau1=-999;
    ptTau2=-999; etaTau2=-999; phiTau2=-999; mTau2=-999;
    ptB1=-999; etaB1=-999; phiB1=-999; mB1=-999;
    ptB2=-999; etaB2=-999; phiB2=-999; mB2=-999;

    ptTau1_gen=-999; etaTau1_gen=-999; phiTau1_gen=-999; mTau1_gen=-999; tauIso1=0;
    ptTau2_gen=-999; etaTau2_gen=-999; phiTau2_gen=-999; mTau2_gen=-999; tauIso2=0;
    ptB1_gen=-999; etaB1_gen=-999; phiB1_gen=-999; mB1_gen=-999;
    ptB2_gen=-999; etaB2_gen=-999; phiB2_gen=-999; mB2_gen=-999;
    ptC1_gen=-999; etaC1_gen=-999; phiC1_gen=-999; mC1_gen=-999;
    ptC2_gen=-999; etaC2_gen=-999; phiC2_gen=-999; mC2_gen=-999;
    ptJ1_gen=-999; etaJ1_gen=-999; phiJ1_gen=-999; mJ1_gen=-999;
    ptJ2_gen=-999; etaJ2_gen=-999; phiJ2_gen=-999; mJ2_gen=-999;
    mTT=-999; mBB1=-999; mt2pileup=-999;
    eventType=2;

    //eventWeight=1.24*341*1000; //fb^-1
    eventWeight=2.92;

    metPx=info->metPx;
    metPy=info->metPy;

    TVector2 genMet(metPx, metPy);

    Float_t sMetPx=rng->Gaus(info->metPx, 20);
    Float_t sMetPy=rng->Gaus(info->metPy, 20);

    TVector2 smearedMet(sMetPx, sMetPy);

    met2=smearedMet.Mod();
    met2Phi=smearedMet.Phi();

    Int_t isLep=0;

    for (Int_t j=0; j<part->GetEntries(); j++) {
      const baconhep::TGenParticle* genloop = (baconhep::TGenParticle*) ((*part)[j]);

      Int_t parentPdg=dynamic_cast<baconhep::TGenParticle *>(part->At(genloop->parent>-1 ? genloop->parent : 0))->pdgId;
      
      if (!(fabs(genloop->pdgId)==15)&&!(fabs(parentPdg)==15)) continue;

      if ( (fabs(genloop->pdgId)==13||fabs(genloop->pdgId)==11) && (fabs(parentPdg)==15) ) {
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

    //if (isLep==1 || iTau1==-1 || iTau2==-1) continue;

    const baconhep::TGenParticle* genTau1 = (baconhep::TGenParticle*) ((*part)[iTau1]);
    const baconhep::TGenParticle* genTau2 = (baconhep::TGenParticle*) ((*part)[iTau2]);
    
    LorentzVector vGenTau1(genTau1->pt,genTau1->eta,genTau1->phi,genTau1->mass);
    LorentzVector vGenTau2(genTau2->pt,genTau2->eta,genTau2->phi,genTau2->mass);

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

    LorentzVector genSumPt(0,0,0,0);
    LorentzVector smeSumPt(0,0,0,0);

    for (Int_t j=0; j<jet->GetEntries(); j++) {
      const baconhep::TGenJet* loop = (baconhep::TGenJet*) ((*jet)[j]);
      LorentzVector vTempG(loop->pt, loop->eta, loop->phi, loop->mass);
      LorentzVector vTempS(loop->pt*(rng->Gaus(1,res(loop->pt, loop->eta))), loop->eta, loop->phi, loop->mass);
      genSumPt+=vTempG;
      smeSumPt+=vTempS;
      
      if (fabs(loop->mcFlavor)==5 && iB1==-1) iB1=j;
      else if (fabs(loop->mcFlavor)==5 && iB2==-1) iB2=j;

      if (iTau1>-1) {
	if (deltaR(loop->eta,vGenTau1.Eta(),loop->phi,vGenTau1.Phi())<0.4) iT1=j;
      }
      if (iTau2>-1) {
	if (deltaR(loop->eta,vGenTau2.Eta(),loop->phi,vGenTau2.Phi())<0.4) iT2=j;
      }
    }
    if (iB1>-1) {
      const baconhep::TGenJet *bjet = (baconhep::TGenJet*) ((*jet)[iB1]);
      ptB1_gen=bjet->pt;
      ptB1=ptB1_gen*(rng->Gaus(1,res(bjet->pt, bjet->eta)));
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
      ptB2=ptB2_gen*(rng->Gaus(1,res(bjet->pt, bjet->eta)));
      etaB2_gen=bjet->eta;
      etaB2=bjet->eta;
      mB2_gen=bjet->mass;
      mB2=bjet->mass;
      phiB2_gen=bjet->phi;
      phiB2=bjet->phi;
      eventWeight*=btag(ptB2, etaB2);
    }

    /*    if (iT1==-1 || iT2==-1) continue;

    const baconhep::TGenJet *taujet = (baconhep::TGenJet*) ((*jet)[iT1]);
    const baconhep::TGenJet *taujet2 = (baconhep::TGenJet*) ((*jet)[iT2]);

    Float_t fact1=(rng->Gaus(1,res(taujet->pt, taujet->eta)));
    Float_t fact2=(rng->Gaus(1,res(taujet2->pt, taujet2->eta)));

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

      ptTau2_gen=taujet->pt;
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
      ptTau2_gen=taujet->pt;
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
    eventWeight*=0.65*0.65;
    */    
    if (ptB1==-999 || ptB2==-999) continue; // || ptTau1==-999 || ptTau2==-999) continue;
    iHH++;
    /*
    LorentzVector metCorr = smeSumPt - genSumPt;
    TVector2 metCorrT(metCorr.Pt(), metCorr.Phi());
    TVector2 newMet=genMet-newMet;
    met=newMet.Mod();
    metPhi=newMet.Phi();

    LorentzVector recoB1(ptB1, etaB1, phiB1, mB1);
    LorentzVector recoB2(ptB2, etaB2, phiB2, mB2);
    LorentzVector recoTau1(ptTau1, etaTau1, phiTau1, mTau1);
    LorentzVector recoTau2(ptTau2, etaTau2, phiTau2, mTau2);
    
    LorentzVector vTT=recoTau1+recoTau2;
    mTT=vTT.M();

    LorentzVector vBB=recoB1+recoB2;
    mBB1=vBB.M();

    //if (recoTau1.Pt()>0 && recoTau2.Pt()>0 && recoB1.Pt()>0 && recoB2.Pt()>0) {
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
    
    if (recoTau1.Pt()>0 && recoTau2.Pt()>0) {
      int channel=0;
      mithep::TSVfit svfit;
      svfit.cov_00=lcov00pp;
      svfit.cov_01=lcov01pp;
      svfit.cov_10=lcov10pp;
      svfit.cov_11=lcov11pp;
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
      m_svpileup = fitter->integrate(&svfit,met2,met2Phi,channel);
      std::cout << tauCat1 << "  " <<  tauCat2 << "   " << channel  << "  "  << m_sv << "  "  <<  m_svpileup << "  " << std::endl; 
    }
    */
    outtree->Fill();
    
  }

  cout << "hh " << iHH << endl;
  cout << "eh " << iEH << endl;
  cout << "mh " << iMH << endl;

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

Double_t res(Double_t pt, Double_t eta) {

  if (fabs(eta)<3.0) return 0.15;
  else if (fabs(eta)<5.0) return 0.3;
  else return 0;

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
