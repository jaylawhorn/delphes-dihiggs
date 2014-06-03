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
#include <THStack.h>
#include <TLine.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <vector>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <TGaxis.h>
#include <TLorentzVector.h>
#include "Math/LorentzVector.h"
#include "yieldData.hh"
#endif

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;

void confParse(const TString conf, vector<TString> &sampleNames, vector<TString> &sampleTitles, vector<Int_t> &sampleColors);

Float_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 );
Double_t deltaPhi(const Double_t phi1, const Double_t phi2);

yieldData getYield(const TString infilename, Float_t TAU_PT_MIN_H, Float_t TAU_PT_MIN_L, Float_t TAU_PT_MIN_HH, Float_t B_PT_MIN, Float_t THRESH, Float_t TAU_ETA_MAX, Float_t B_ETA_MAX);

void cat(Float_t THRESH=50) {

  Float_t TAU_PT_MIN_H = 30;
  Float_t TAU_PT_MIN_L = 30;
  Float_t TAU_PT_MIN_HH = 30;
  Float_t B_PT_MIN = 30;
  
  Float_t TAU_ETA_MAX = 4.0;
  Float_t B_ETA_MAX = 4.0;

  yieldData signal = getYield("/afs/cern.ch/work/j/jlawhorn/public/ntuples/HHToTTBB_14TeV.root", TAU_PT_MIN_H, TAU_PT_MIN_L, TAU_PT_MIN_HH, B_PT_MIN, THRESH, TAU_ETA_MAX, B_ETA_MAX);
  yieldData ttbar = getYield("/afs/cern.ch/work/j/jlawhorn/public/ntuples/tt.root", TAU_PT_MIN_H, TAU_PT_MIN_L, TAU_PT_MIN_HH, B_PT_MIN, THRESH, TAU_ETA_MAX, B_ETA_MAX);
  yieldData bkgd = getYield("/afs/cern.ch/work/j/jlawhorn/public/ntuples/bkgd.root", TAU_PT_MIN_H, TAU_PT_MIN_L, TAU_PT_MIN_HH, B_PT_MIN, THRESH, TAU_ETA_MAX, B_ETA_MAX);

  cout << setprecision(1) << fixed;

  cout << "                ALL      2BTAG      1BTAG      0BTAG " << endl;
  cout << "tt (sig) " << setw(10) << signal.total_tt*3e3 << " " << setw(10) << signal.tt_2*3e3 << " " << setw(10) <<  signal.tt_1*3e3 << " " << setw(10) << signal.tt_0*3e3 << endl; 
  cout << " (ttbar) " << setw(10) << ttbar.total_tt*3e3 << " " << setw(10) << ttbar.tt_2*3e3 << " " << setw(10) << ttbar.tt_1*3e3 << " " << setw(10) << ttbar.tt_0*3e3 << endl; 
  cout << " (other) " << setw(10) << bkgd.total_tt*3e3 << " " << setw(10) << bkgd.tt_2*3e3 << " " << setw(10) << bkgd.tt_1*3e3 << " " << setw(10) << bkgd.tt_0*3e3 << endl; 
  cout << " ----" << endl;
  cout << "mt (sig) " << setw(10) << signal.total_mt*3e3 << " " << setw(10) << signal.mt_2*3e3 << " " << setw(10) << signal.mt_1*3e3 << " " << setw(10) << signal.mt_0*3e3 << endl; 
  cout << " (ttbar) " << setw(10) << ttbar.total_mt*3e3 << " " << setw(10) << ttbar.mt_2*3e3 << " " << setw(10) << ttbar.mt_1*3e3 << " " << setw(10) << ttbar.mt_0*3e3 << endl; 
  cout << " (other) " << setw(10) << bkgd.total_mt*3e3 << " " << setw(10) << bkgd.mt_2*3e3 << " " << setw(10) << bkgd.mt_1*3e3 << " " << setw(10) << bkgd.mt_0*3e3 << endl; 
  cout << " ----" << endl;
  cout << "et (sig) " << setw(10) << signal.total_et*3e3 << " " << setw(10) << signal.et_2*3e3 << " " << setw(10) << signal.et_1*3e3 << " " << setw(10) << signal.et_0*3e3 << endl; 
  cout << " (ttbar) " << setw(10) << ttbar.total_et*3e3 << " " << setw(10) << ttbar.et_2*3e3 << " " << setw(10) << ttbar.et_1*3e3 << " " << setw(10) << ttbar.et_0*3e3 << endl; 
  cout << " (other) " << setw(10) << bkgd.total_et*3e3 << " " << setw(10) << bkgd.et_2*3e3 << " " << setw(10) << bkgd.et_1*3e3 << " " << setw(10) << bkgd.et_0*3e3 << endl; 
  cout << " ----" << endl;
  cout << "em (sig) " << setw(10) << signal.total_em*3e3 << " " << setw(10) << signal.em_2*3e3 << " " << setw(10) << signal.em_1*3e3 << " " << setw(10) << signal.em_0*3e3 << endl; 
  cout << " (ttbar) " << setw(10) << ttbar.total_em*3e3 << " " << setw(10) << ttbar.em_2*3e3 << " " << setw(10) << ttbar.em_1*3e3 << " " << setw(10) << ttbar.em_0*3e3 << endl; 
  cout << " (other) " << setw(10) << bkgd.total_em*3e3 << " " << setw(10) << bkgd.em_2*3e3 << " " << setw(10) << bkgd.em_1*3e3 << " " << setw(10) << bkgd.em_0*3e3 << endl; 
  cout << endl;
  cout << " BOOSTED REQ " << endl;
  cout << "              2BTAG      1BTAG      0BTAG " << endl;
  cout << "tt (sig) " << setw(10) << signal.tt_2_boo*3e3 << " " << setw(10) <<  signal.tt_1_boo*3e3 << " " << setw(10) << signal.tt_0_boo*3e3 << endl; 
  cout << " (ttbar) " << setw(10) << ttbar.tt_2_boo*3e3 << " " << setw(10) << ttbar.tt_1_boo*3e3 << " " << setw(10) << ttbar.tt_0_boo*3e3 << endl; 
  cout << " (other) " << setw(10) << setw(10) << bkgd.tt_2_boo*3e3 << " " << setw(10) << bkgd.tt_1_boo*3e3 << " " << setw(10) << bkgd.tt_0_boo*3e3 << endl; 
  cout << " ----" << endl;
  cout << "mt (sig) " << setw(10) << signal.mt_2_boo*3e3 << " " << setw(10) << signal.mt_1_boo*3e3 << " " << setw(10) << signal.mt_0_boo*3e3 << endl; 
  cout << " (ttbar) " << setw(10) << ttbar.mt_2_boo*3e3 << " " << setw(10) << ttbar.mt_1_boo*3e3 << " " << setw(10) << ttbar.mt_0_boo*3e3 << endl; 
  cout << " (other) " << setw(10) << bkgd.mt_2_boo*3e3 << " " << setw(10) << bkgd.mt_1_boo*3e3 << " " << setw(10) << bkgd.mt_0_boo*3e3 << endl; 
  cout << " ----" << endl;
  cout << "et (sig) " << setw(10) << signal.et_2_boo*3e3 << " " << setw(10) << signal.et_1_boo*3e3 << " " << setw(10) << signal.et_0_boo*3e3 << endl; 
  cout << " (ttbar) " << setw(10) << ttbar.et_2_boo*3e3 << " " << setw(10) << ttbar.et_1_boo*3e3 << " " << setw(10) << ttbar.et_0_boo*3e3 << endl; 
  cout << " (other) " << setw(10) << bkgd.et_2_boo*3e3 << " " << setw(10) << bkgd.et_1_boo*3e3 << " " << setw(10) << bkgd.et_0_boo*3e3 << endl; 
  cout << " ----" << endl;
  cout << "em (sig) " << setw(10) << signal.em_2_boo*3e3 << " " << setw(10) << signal.em_1_boo*3e3 << " " << setw(10) << signal.em_0_boo*3e3 << endl; 
  cout << " (ttbar) " << setw(10) << ttbar.em_2_boo*3e3 << " " << setw(10) << ttbar.em_1_boo*3e3 << " " << setw(10) << ttbar.em_0_boo*3e3 << endl; 
  cout << " (other) " << setw(10) << bkgd.em_2_boo*3e3 << " " << setw(10) << bkgd.em_1_boo*3e3 << " " << setw(10) << bkgd.em_0_boo*3e3 << endl; 
  cout << endl;
  cout << " NOT BOOSTED " << endl;
  cout << "              2BTAG      1BTAG      0BTAG " << endl;
  cout << "tt (sig) " << setw(10) << signal.tt_2_nb*3e3 << " " << setw(10) <<  signal.tt_1_nb*3e3 << " " << setw(10) << signal.tt_0_nb*3e3 << endl; 
  cout << " (ttbar) " << setw(10) << ttbar.tt_2_nb*3e3 << " " << setw(10) << ttbar.tt_1_nb*3e3 << " " << setw(10) << ttbar.tt_0_nb*3e3 << endl; 
  cout << " (other) " << setw(10) << setw(10) << bkgd.tt_2_nb*3e3 << " " << setw(10) << bkgd.tt_1_nb*3e3 << " " << setw(10) << bkgd.tt_0_nb*3e3 << endl; 
  cout << " ----" << endl;
  cout << "mt (sig) " << setw(10) << signal.mt_2_nb*3e3 << " " << setw(10) << signal.mt_1_nb*3e3 << " " << setw(10) << signal.mt_0_nb*3e3 << endl; 
  cout << " (ttbar) " << setw(10) << ttbar.mt_2_nb*3e3 << " " << setw(10) << ttbar.mt_1_nb*3e3 << " " << setw(10) << ttbar.mt_0_nb*3e3 << endl; 
  cout << " (other) " << setw(10) << bkgd.mt_2_nb*3e3 << " " << setw(10) << bkgd.mt_1_nb*3e3 << " " << setw(10) << bkgd.mt_0_nb*3e3 << endl; 
  cout << " ----" << endl;
  cout << "et (sig) " << setw(10) << signal.et_2_nb*3e3 << " " << setw(10) << signal.et_1_nb*3e3 << " " << setw(10) << signal.et_0_nb*3e3 << endl; 
  cout << " (ttbar) " << setw(10) << ttbar.et_2_nb*3e3 << " " << setw(10) << ttbar.et_1_nb*3e3 << " " << setw(10) << ttbar.et_0_nb*3e3 << endl; 
  cout << " (other) " << setw(10) << bkgd.et_2_nb*3e3 << " " << setw(10) << bkgd.et_1_nb*3e3 << " " << setw(10) << bkgd.et_0_nb*3e3 << endl; 
  cout << " ----" << endl;
  cout << "em (sig) " << setw(10) << signal.em_2_nb*3e3 << " " << setw(10) << signal.em_1_nb*3e3 << " " << setw(10) << signal.em_0_nb*3e3 << endl; 
  cout << " (ttbar) " << setw(10) << ttbar.em_2_nb*3e3 << " " << setw(10) << ttbar.em_1_nb*3e3 << " " << setw(10) << ttbar.em_0_nb*3e3 << endl; 
  cout << " (other) " << setw(10) << bkgd.em_2_nb*3e3 << " " << setw(10) << bkgd.em_1_nb*3e3 << " " << setw(10) << bkgd.em_0_nb*3e3 << endl; 

}

yieldData getYield(const TString infilename, Float_t TAU_PT_MIN_H, Float_t TAU_PT_MIN_L, Float_t TAU_PT_MIN_HH, Float_t B_PT_MIN, Float_t THRESH, Float_t TAU_ETA_MAX, Float_t B_ETA_MAX) {

  // tau decay modes
  enum { hadron=1, electron, muon };

  // tautau decay modes
  enum { all=0, dijet, jetmu, jetele, muele };

  yieldData data;

  data.total_yield=0;
  data.total_tt=0; data.total_mt=0; data.total_et=0; data.total_em=0;
  data.tt_2=0; data.mt_2=0; data.et_2=0; data.em_2=0; 
  data.tt_1=0; data.mt_1=0; data.et_1=0; data.em_1=0; 
  data.tt_0=0; data.mt_0=0; data.et_0=0; data.em_0=0; 

  data.tt_2_boo=0; data.mt_2_boo=0; data.et_2_boo=0; data.em_2_boo=0; 
  data.tt_1_boo=0; data.mt_1_boo=0; data.et_1_boo=0; data.em_1_boo=0; 
  data.tt_0_boo=0; data.mt_0_boo=0; data.et_0_boo=0; data.em_0_boo=0; 

  data.tt_2_nb=0; data.mt_2_nb=0; data.et_2_nb=0; data.em_2_nb=0; 
  data.tt_1_nb=0; data.mt_1_nb=0; data.et_1_nb=0; data.em_1_nb=0; 
  data.tt_0_nb=0; data.mt_0_nb=0; data.et_0_nb=0; data.em_0_nb=0; 

  Float_t eventWeight=1;
  UInt_t eventType;
  UInt_t bTag1, bTag2;
  UInt_t tauCat1, tauCat2;
  UInt_t tFake1, tFake2;
  UInt_t bFake1, bFake2;
  Float_t met, metPhi;
  Double_t mt2;
  LorentzVector *sRecoB1=0, *sRecoB2=0;
  LorentzVector *sGenB1=0, *sGenB2=0;
  LorentzVector *sGenJetB1=0, *sGenJetB2=0;
  LorentzVector *sRecoTau1=0, *sRecoTau2=0;
  LorentzVector *sGenTau1=0, *sGenTau2=0;
  LorentzVector *sGenJetTau1=0, *sGenJetTau2=0;
  LorentzVector *sRecoJet=0;

  TFile *infile;
  TTree *intree;

  //cout << "Processing  " << infilename << " ..." << endl;
  infile = new TFile(infilename); assert(infile);
  intree = (TTree*) infile->Get("Events"); assert(intree);
  
  intree->SetBranchAddress("eventWeight",    &eventWeight);
  intree->SetBranchAddress("eventType",      &eventType);
  intree->SetBranchAddress("tauCat1",        &tauCat1);
  intree->SetBranchAddress("tauCat2",        &tauCat2);
  intree->SetBranchAddress("bTag1",          &bTag1);
  intree->SetBranchAddress("bTag2",          &bTag2);
  intree->SetBranchAddress("tFake1",         &tFake1);
  intree->SetBranchAddress("tFake2",         &tFake2);
  intree->SetBranchAddress("bFake1",         &bFake1);
  intree->SetBranchAddress("bFake2",         &bFake2);
  intree->SetBranchAddress("met",            &met);
  intree->SetBranchAddress("metPhi",         &metPhi);
  intree->SetBranchAddress("mt2",            &mt2);
  intree->SetBranchAddress("sGenTau1",       &sGenTau1);      // 4-vector for generator leading tau
  intree->SetBranchAddress("sGenTau2",       &sGenTau2);      // 4-vector for generator second tau
  intree->SetBranchAddress("sGenB1",         &sGenB1);        // 4-vector for generator leading b-jet
  intree->SetBranchAddress("sGenB2",         &sGenB2);        // 4-vector for generator second b-jet
  intree->SetBranchAddress("sGenJetTau1",    &sGenJetTau1);   // 4-vector for generator leading tau
  intree->SetBranchAddress("sGenJetTau2",    &sGenJetTau2);   // 4-vector for generator second tau
  intree->SetBranchAddress("sGenJetB1",      &sGenJetB1);     // 4-vector for generator leading b-jet
  intree->SetBranchAddress("sGenJetB2",      &sGenJetB2);     // 4-vector for generator second b-jet
  intree->SetBranchAddress("sRecoTau1",      &sRecoTau1);     // 4-vector for reconstructed leading tau
  intree->SetBranchAddress("sRecoTau2",      &sRecoTau2);     // 4-vector for reconstructed second tau
  intree->SetBranchAddress("sRecoB1",        &sRecoB1);       // 4-vector for reconstructed leading b-jet
  intree->SetBranchAddress("sRecoB2",        &sRecoB2);       // 4-vector for reconstructed second b-jet
  intree->SetBranchAddress("sRecoJet",       &sRecoJet);      // 4-vector for reconstructed extra jet

  Int_t CAT=-1, BOOST=-1, NB=0;
  
  for(UInt_t iEntry=0; iEntry<intree->GetEntries(); iEntry++) { // entry loop
    intree->GetEntry(iEntry);

    //eventWeight=1;

    CAT=-1; BOOST=-1; NB=0;

    // skip events that don't have 2 reco b's and 2 reco tau's
    if ( (sRecoB1->Pt()==999) || (sRecoB2->Pt()==999) ) continue;
    if ( (sRecoTau1->Pt()==999) || (sRecoTau2->Pt()==999) ) continue;
    
    data.total_yield+=eventWeight;

    if ((sRecoB1->Pt()<B_PT_MIN)||(sRecoB2->Pt()<B_PT_MIN)) continue;
    if ((fabs(sRecoB1->Eta())>B_ETA_MAX)||(fabs(sRecoB2->Eta())>B_ETA_MAX)) continue;

    if ((fabs(sRecoTau1->Eta())>TAU_ETA_MAX)||(fabs(sRecoTau2->Eta())>TAU_ETA_MAX)) continue;

    if ( ( (tauCat1==electron) && (tauCat2==electron) ) || ( (tauCat1==muon) && (tauCat2==muon) ) ) continue;
    else if ( ( tauCat1 == hadron ) && ( tauCat2==hadron ) ) {
      if ( (sRecoTau1->Pt()<TAU_PT_MIN_H)||(sRecoTau2->Pt()<TAU_PT_MIN_HH) ) continue;
      CAT=2;
    }
    else if ( ( tauCat1==hadron ) && ( tauCat2==muon) ) {
      if ( (sRecoTau1->Pt()<TAU_PT_MIN_H)||(sRecoTau2->Pt()<TAU_PT_MIN_L) ) continue;
      CAT=3;
    }
    else if ( ( tauCat2==hadron ) && ( tauCat1==muon) ) {
      if ( (sRecoTau1->Pt()<TAU_PT_MIN_L)||(sRecoTau2->Pt()<TAU_PT_MIN_H) ) continue;
      CAT=3;
    }
    else if ( (tauCat1==hadron) && (tauCat2==electron) ) {
      if ( (sRecoTau1->Pt()<TAU_PT_MIN_H)||(sRecoTau2->Pt()<TAU_PT_MIN_L) ) continue;
      CAT=4;
    }
    else if ( (tauCat2==hadron) && (tauCat1==electron) ) {
      if ( (sRecoTau1->Pt()<TAU_PT_MIN_L)||(sRecoTau2->Pt()<TAU_PT_MIN_H) ) continue;
      CAT=4;
    }
    else if ( ( (tauCat1==muon) && (tauCat2==electron) ) || ( (tauCat2==muon) && (tauCat1==electron) ) ) {
      if ( (sRecoTau1->Pt()<TAU_PT_MIN_L)||(sRecoTau2->Pt()<TAU_PT_MIN_L) ) continue;
      CAT=5;
    }

    if (bTag1==3) bTag1=2;
    if (bTag2==3) bTag2=2;

    if ( (bTag1==2) && (bTag2==2)) NB=2;
    else if ( (bTag1>0) || (bTag2>0) ) NB=1;
    else if ( (bTag1==0) || (bTag2==0) ) NB=0;
    else cout << "btagging is messed up?" << endl;

    //if ((sRecoJet->Pt()<THRESH)||(sRecoJet->Pt()==999)) {
    //BOOST=0;
    //}
    //else if (sRecoJet->Pt()>THRESH) BOOST=1;
    //if (THRESH<1) BOOST=1;

    LorentzVector hTT=*sRecoTau1+*sRecoTau2;
    LorentzVector hBB=*sRecoB1+*sRecoB2;

    if ((hTT.M()<50)||(hTT.M()>150)) continue;
    if ((hBB.M()<100)||(hBB.M()>150)) continue;
    
    LorentzVector hh=hTT+hBB;

    if (hh.M()<300) continue;    

    if (mt2<150) continue;

    if (hh.Pt()<THRESH) BOOST=0;
    else if (hh.Pt()>THRESH) BOOST=1;
    
    if (THRESH<1) BOOST=1;

    // TAU-TAU
    if (CAT==2) {
      data.total_tt+=eventWeight; 

      if (NB==2) {
	data.tt_2+=eventWeight;
	if (BOOST==1) data.tt_2_boo+=eventWeight;
	else data.tt_2_nb+=eventWeight;
      }
      else if (NB==1) {
	data.tt_1+=eventWeight;
	if (BOOST==1) data.tt_1_boo+=eventWeight;
	else data.tt_1_nb+=eventWeight;
      }
      else if (NB==0) {
	data.tt_0+=eventWeight;
	if (BOOST==1) data.tt_0_boo+=eventWeight;
	else data.tt_0_nb+=eventWeight;
      }
    }

    // MU-TAU
    else if (CAT==3) {
      data.total_mt+=eventWeight; 

      if (NB==2) {
	data.mt_2+=eventWeight;
	if (BOOST==1) data.mt_2_boo+=eventWeight;
	else data.mt_2_nb+=eventWeight;
      }
      else if (NB==1) {
	data.mt_1+=eventWeight;
	if (BOOST==1) data.mt_1_boo+=eventWeight;
	else data.mt_1_nb+=eventWeight;
      }
      else if (NB==0) {
	data.mt_0+=eventWeight;
	if (BOOST==1) data.mt_0_boo+=eventWeight;
	else data.mt_0_nb+=eventWeight;
      }
    }

    // ELE-TAU
    else if (CAT==4) {
      data.total_et+=eventWeight; 

      if (NB==2) {
	data.et_2+=eventWeight;
	if (BOOST==1) data.et_2_boo+=eventWeight;
	else data.et_2_nb+=eventWeight;
      }
      else if (NB==1) {
	data.et_1+=eventWeight;
	if (BOOST==1) data.et_1_boo+=eventWeight;
	else data.et_1_nb+=eventWeight;
      }
      else if (NB==0) {
	data.et_0+=eventWeight;
	if (BOOST==1) data.et_0_boo+=eventWeight;
	else data.et_0_nb+=eventWeight;
      }
    }

    // ELE-MU
    else if (CAT==5) {
      data.total_em+=eventWeight; 

      if (NB==2) {
	data.em_2+=eventWeight;
	if (BOOST==1) data.em_2_boo+=eventWeight;
	else data.em_2_nb+=eventWeight;
      }
      else if (NB==1) {
	data.em_1+=eventWeight;
	if (BOOST==1) data.em_1_boo+=eventWeight;
	else data.em_1_nb+=eventWeight;
      }
      else if (NB==0) {
	data.em_0+=eventWeight;
	if (BOOST==1) data.em_0_boo+=eventWeight;
	else data.em_0_nb+=eventWeight;
      }
    }

  } // end entry loop
  
  delete infile;
  infile=0, intree=0;

  return data;

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

Double_t deltaPhi(Double_t phi1, Double_t phi2) 
{
  // Compute dPhi between two given angles. Results is in [0,pi].
  const Double_t pi = 3.14159265358979;
  Double_t dphi = fabs(phi1-phi2);
  while (dphi>pi)
    dphi = fabs(dphi - 2.0*pi);

  return dphi;
}

Float_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 ) {

  const Float_t pi = 3.14159265358979;

  Float_t etaDiff = (eta1-eta2);
  Float_t phiDiff = fabs(phi1-phi2);
  while ( phiDiff>pi ) phiDiff = fabs(phiDiff-2.0*pi);

  Float_t deltaRSquared = etaDiff*etaDiff + phiDiff*phiDiff;

  return TMath::Sqrt(deltaRSquared);

}
