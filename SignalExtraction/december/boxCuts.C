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
#include "../../MitStyleRemix.hh"
#include "../CPlot.hh"

#endif

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;

void confParse(const TString conf, vector<TString> &sampleNames, vector<TString> &sampleTitles, vector<Int_t> &sampleColors);

Float_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 );
Double_t deltaPhi(const Double_t phi1, const Double_t phi2);

void boxCuts(const TString conf="new.conf", Int_t tauType=0) {

  cout << "Selected ditau decay: ";
  if (tauType==0) { 
    cout << "all" << endl;
    CPlot::sOutDir = "all";     
  }
  else if (tauType==1) {
    cout << "dijet" << endl;
    CPlot::sOutDir = "dijet"; 
  }
  else if (tauType==2) {
    cout << "jet-mu" << endl;
    CPlot::sOutDir = "jetmu"; 
  }
  else if (tauType==3) {
    cout << "jet-ele" << endl;
      CPlot::sOutDir = "jetele"; 
  }
  else if (tauType==4) {
    cout << "mu-ele" << endl;
    CPlot::sOutDir = "muele"; 
  }

  // define pre-selection variables
  const Float_t TAU_PT_MIN_H = 30;
  const Float_t TAU_PT_MIN_L = 20;
  const Float_t B_PT_MIN = 30;
  
  const Float_t TAU_ETA_MAX = 2.5;
  const Float_t B_ETA_MAX = 2.5;
  const Int_t ETA_BINS = 12;

  // tau decay modes
  enum { hadron=1, electron, muon };

  vector<TString> sampleNames;
  vector<TString> sampleTitles;
  vector<Int_t> sampleColors;

  Double_t hhTotal=0, ttTotal=0, bbTotal=0, llTotal=0;
  Double_t hhPre=0, ttPre=0, bbPre=0, llPre=0;
  Double_t hh1=0, tt1=0, bb1=0, ll1=0;
  Double_t hh2=0, tt2=0, bb2=0, ll2=0;
  Double_t hh3=0, tt3=0, bb3=0, ll3=0;

  Double_t hhCount=0, ttCount=0, bbCount=0, llCount=0;

  confParse(conf, sampleNames, sampleTitles, sampleColors);

  vector<TH1D*> hTpt1, hTpt2, hBpt1, hBpt2;
  vector<TH1D*> hTeta1, hTeta2, hBeta1, hBeta2;
  vector<TH1D*> hTTMass, hBBMass, hTTBBMass;
  vector<TH1D*> hTTAngle, hBBAngle;
  vector<TH1D*> hMt2;

  char hname[100];
  
  for(UInt_t isam=0; isam<sampleNames.size(); isam++) {
    sprintf(hname, "hTpt1_%i",isam); hTpt1.push_back(new TH1D(hname, sampleTitles[isam], 100, 0, 300)); hTpt1[isam]->Sumw2(); 
    sprintf(hname, "hTpt2_%i",isam); hTpt2.push_back(new TH1D(hname, sampleTitles[isam], 100, 0, 300)); hTpt2[isam]->Sumw2(); 
    sprintf(hname, "hBpt1_%i",isam); hBpt1.push_back(new TH1D(hname, sampleTitles[isam], 100, 0, 300)); hBpt1[isam]->Sumw2(); 
    sprintf(hname, "hBpt2_%i",isam); hBpt2.push_back(new TH1D(hname, sampleTitles[isam], 100, 0, 300)); hBpt2[isam]->Sumw2(); 

    sprintf(hname, "hTeta1_%i",isam); hTeta1.push_back(new TH1D(hname, sampleTitles[isam], 16, -2.5, 2.5)); hTeta1[isam]->Sumw2(); 
    sprintf(hname, "hTeta2_%i",isam); hTeta2.push_back(new TH1D(hname, sampleTitles[isam], 16, -2.5, 2.5)); hTeta2[isam]->Sumw2(); 
    sprintf(hname, "hBeta1_%i",isam); hBeta1.push_back(new TH1D(hname, sampleTitles[isam], 16, -2.5, 2.5)); hBeta1[isam]->Sumw2(); 
    sprintf(hname, "hBeta2_%i",isam); hBeta2.push_back(new TH1D(hname, sampleTitles[isam], 16, -2.5, 2.5)); hBeta2[isam]->Sumw2(); 

    sprintf(hname, "hTTMass_%i",isam); hTTMass.push_back(new TH1D(hname, sampleTitles[isam], 100, 0, 400)); hTTMass[isam]->Sumw2(); 
    sprintf(hname, "hBBMass_%i",isam); hBBMass.push_back(new TH1D(hname, sampleTitles[isam], 100, 0, 400)); hBBMass[isam]->Sumw2(); 
    sprintf(hname, "hTTBBMass_%i",isam); hTTBBMass.push_back(new TH1D(hname, sampleTitles[isam], 100, 0, 1200)); hTTBBMass[isam]->Sumw2(); 
    sprintf(hname, "hTTAngle_%i",isam); hTTAngle.push_back(new TH1D(hname, sampleTitles[isam], 20, 0, TMath::Pi())); hTTAngle[isam]->Sumw2(); 
    sprintf(hname, "hBBAngle_%i",isam); hBBAngle.push_back(new TH1D(hname, sampleTitles[isam], 20, 0, TMath::Pi())); hBBAngle[isam]->Sumw2(); 
    sprintf(hname, "hMt2_%i",isam); hMt2.push_back(new TH1D(hname, sampleTitles[isam], 50, 0, 500)); hMt2[isam]->Sumw2(); 
  }

  Float_t eventWeight=1;
  UInt_t bTag1, bTag2;
  UInt_t tauCat1, tauCat2;
  Float_t met, metPhi;
  Double_t mt2;
  LorentzVector *recoB1=0, *recoB2=0;
  LorentzVector *recoTau1=0, *recoTau2=0;
  LorentzVector *recoExtraJet=0;
  LorentzVector *tauHiggs=0, *bHiggs=0;

  TFile *infile;
  TTree *intree;

  Float_t bPt1=0, bPt2=0;
  Float_t tauPt1=0, tauPt2=0;
  Float_t tauEta1=0, tauEta2=0;

  for (UInt_t isam=0; isam<sampleNames.size(); isam++) { // sample loop

    TString infilename = sampleNames[isam];
    //cout << "Processing  " << infilename << " ..." << endl;
    infile = new TFile(infilename); assert(infile);
    intree = (TTree*) infile->Get("Events"); assert(intree);

    intree->SetBranchAddress("eventWeight",    &eventWeight);
    intree->SetBranchAddress("tauCat1",        &tauCat1);
    intree->SetBranchAddress("tauCat2",        &tauCat2);
    intree->SetBranchAddress("bTag1",          &bTag1);
    intree->SetBranchAddress("bTag2",          &bTag2);
    intree->SetBranchAddress("met",            &met);
    intree->SetBranchAddress("metPhi",         &metPhi);
    intree->SetBranchAddress("mt2",            &mt2);
    intree->SetBranchAddress("recoTau1",       &recoTau1);     // 4-vector for reconstructed leading tau
    intree->SetBranchAddress("recoTau2",       &recoTau2);     // 4-vector for reconstructed second tau
    intree->SetBranchAddress("recoB1",         &recoB1);       // 4-vector for reconstructed leading b-jet
    intree->SetBranchAddress("recoB2",         &recoB2);       // 4-vector for reconstructed second b-jet
    intree->SetBranchAddress("recoExtraJet",   &recoExtraJet); // 4-vector for reconstructed extra jet

    for(UInt_t iEntry=0; iEntry<intree->GetEntries(); iEntry++) { // entry loop
      intree->GetEntry(iEntry);

      if (isam==0) hhTotal+=eventWeight;
      if (isam==1) ttTotal+=eventWeight;
      if (isam==2) bbTotal+=eventWeight;
      if (isam==3) llTotal+=eventWeight;

      // skip events that don't have 2 reco b's and 2 reco tau's
      if ( (recoB1->Pt()==999) || (recoB2->Pt()==999) ) continue;
      if ( (recoTau1->Pt()==999) || (recoTau2->Pt()==999) ) continue;

      // require boost of at least 50 GeV
      if ( (recoExtraJet->Pt()!=999) && (recoExtraJet->Pt()>50) ) continue;

      // exclude ditau->ee/mumu channels
      if ( (tauCat1 == electron) && (tauCat2 == electron) ) continue;
      if ( (tauCat1 == muon) && (tauCat2 == muon) ) continue;

      // select tautypes

      // dijet                                                                                                                                         
      if ( (tauType==1) && ( ( tauCat1 != hadron ) || ( tauCat2!=hadron ) ) ) continue;

      // jetmu                                                                                                                                         
      if ( (tauType==2) && !( ( ( tauCat1==hadron ) && ( tauCat2==muon) ) || ( ( tauCat2==hadron ) && ( tauCat1==muon) ) ) ) continue;

      // jetele                                                                                                                                        
      if ( (tauType==3) && !( ( (tauCat1==hadron) && (tauCat2==electron) ) || ( (tauCat2==hadron) && (tauCat1==electron) ) ) ) continue;

      //muele                                                                                                                                          
      if ( (tauType==4) && !( ( (tauCat1==muon) && (tauCat2==electron) ) || ( (tauCat2==muon) && (tauCat1==electron) ) ) ) continue;

      bPt1=recoB1->Pt(); 
      bPt2=recoB2->Pt();
      if (tauType==2) {
	if (tauCat1==hadron) {
	  tauPt1 = recoTau1->Pt();
	  tauPt2 = recoTau2->Pt();
	  tauEta1 = recoTau1->Eta();
	  tauEta2 = recoTau2->Eta();
	}
	else {
	  tauPt2 = recoTau1->Pt();
	  tauPt1 = recoTau2->Pt();
	  tauEta2 = recoTau1->Eta();
	  tauEta1 = recoTau2->Eta();
	}
      }
      if (tauType==3) {
	if (tauCat1==hadron) {
	  tauPt1 = recoTau1->Pt();
	  tauPt2 = recoTau2->Pt();
	  tauEta1 = recoTau1->Eta();
	  tauEta2 = recoTau2->Eta();
	}
	else {
	  tauPt2 = recoTau1->Pt();
	  tauPt1 = recoTau2->Pt();
	  tauEta2 = recoTau1->Eta();
	  tauEta1 = recoTau2->Eta();
	}
      }
      if (tauType==4) {
	if (tauCat1==muon) {
	  tauPt1 = recoTau1->Pt();
	  tauPt2 = recoTau2->Pt();
	  tauEta1 = recoTau1->Eta();
	  tauEta2 = recoTau2->Eta();
	}
	else {
	  tauPt2 = recoTau1->Pt();
	  tauPt1 = recoTau2->Pt();
	  tauEta2 = recoTau1->Eta();
	  tauEta1 = recoTau2->Eta();
	}
      }
      else {
	tauPt1 = recoTau1->Pt();
	tauPt2 = recoTau2->Pt();
	tauEta1 = recoTau1->Eta();
	tauEta2 = recoTau2->Eta();
      }

      // tau pre-selection
      if ( ( fabs( recoTau1->Eta() ) > TAU_ETA_MAX ) || ( fabs(recoTau2->Eta() ) > TAU_ETA_MAX ) ) continue;

      if (( tauCat1 == hadron ) && (tauPt1 < TAU_PT_MIN_H )) continue;
      else if (( tauCat1 != hadron ) && (tauPt1 < TAU_PT_MIN_L )) continue;

      if (( tauCat2 == hadron ) && (tauPt2 < TAU_PT_MIN_H )) continue;
      else if (( tauCat2 != hadron ) && (tauPt2 < TAU_PT_MIN_L )) continue;

      // b pre-selection
      if ( ( fabs( recoB1->Eta() ) > B_ETA_MAX ) || ( fabs( recoB2->Eta() ) > B_ETA_MAX ) ) continue;

      if ( ( bPt1 < B_PT_MIN ) || ( bPt2 < B_PT_MIN ) ) continue;

      hTpt1[isam]->Fill(tauPt1, eventWeight);
      hTpt2[isam]->Fill(tauPt2, eventWeight);
      hBpt1[isam]->Fill(bPt1, eventWeight);
      hBpt2[isam]->Fill(bPt2, eventWeight);

      hTeta1[isam]->Fill(recoTau1->Eta(), eventWeight);
      hTeta2[isam]->Fill(recoTau2->Eta(), eventWeight);
      hBeta1[isam]->Fill(recoB1->Eta(), eventWeight);
      hBeta2[isam]->Fill(recoB2->Eta(), eventWeight);

      if (isam==0) hhPre+=eventWeight;
      if (isam==1) ttPre+=eventWeight;
      if (isam==2) bbPre+=eventWeight;
      if (isam==3) llPre+=eventWeight;

      LorentzVector vTau1(tauPt1, recoTau1->Eta(), recoTau1->Phi(), recoTau1->M());
      LorentzVector vTau2(tauPt2, recoTau2->Eta(), recoTau2->Phi(), recoTau2->M());
      LorentzVector vTauHiggs = vTau1 + vTau2;
      
      LorentzVector vB1(bPt1, recoB1->Eta(), recoB1->Phi(), recoB1->M());
      LorentzVector vB2(bPt2, recoB2->Eta(), recoB2->Phi(), recoB2->M());
      LorentzVector vBHiggs = vB1 + vB2;

      LorentzVector vHH = vTauHiggs + vBHiggs;

      for (Int_t c=0; c<7; c++) {

	if ( (c!=0) && ((vTauHiggs.M()<50) || (vTauHiggs.M()>120)) ) continue;
	if ( (c!=1) && ((vBHiggs.M()<100) || (vBHiggs.M()>150)) ) continue;
	if ( (c!=2) && (vHH.M()<360) ) continue;
	if ( (c!=3) && (mt2<120) ) continue;
	//if ( (c!=4) && (deltaPhi(recoTau1->Phi(), recoTau2->Phi())<0) ) continue;
	//if ( (c!=4) && (deltaPhi(recoTau1->Phi(), recoTau2->Phi())>1.5) ) continue;
	//if ( (c!=5) && (deltaPhi(recoB1->Phi(), recoB2->Phi())>1.3) ) continue;
	//if ( (c!=5) && (deltaPhi(recoB1->Phi(), recoB2->Phi())<0) ) continue;

	if (c==0) hTTMass[isam]->Fill(vTauHiggs.M(),eventWeight);
	else if (c==1) hBBMass[isam]->Fill(vBHiggs.M(),eventWeight);
	else if (c==2) hTTBBMass[isam]->Fill(vHH.M(),eventWeight);
	else if (c==3) hMt2[isam]->Fill(mt2,eventWeight);
	else if (c==4) hTTAngle[isam]->Fill(deltaPhi(recoTau1->Phi(), recoTau2->Phi()),eventWeight);
	else if (c==5) hBBAngle[isam]->Fill(deltaPhi(recoB1->Phi(), recoB2->Phi()),eventWeight);

      }

      if ( (vTauHiggs.M()<50) || (vTauHiggs.M()>120) ) continue;

      if (isam==0) hh1+=eventWeight;
      if (isam==1) tt1+=eventWeight;
      if (isam==2) bb1+=eventWeight;
      if (isam==3) ll1+=eventWeight;
      
      if ( (vBHiggs.M()<100) || (vBHiggs.M()>150) ) continue;
      
      if (isam==0) hh2+=eventWeight;
      if (isam==1) tt2+=eventWeight;
      if (isam==2) bb2+=eventWeight;
      if (isam==3) ll2+=eventWeight;
      
      if ( vHH.M()<360 ) continue;
      
      if (isam==0) hh3+=eventWeight;
      if (isam==1) tt3+=eventWeight;
      if (isam==2) bb3+=eventWeight;
      if (isam==3) ll3+=eventWeight;
      
      if ( mt2<120 ) continue;
      
      if (isam==0) hhCount+=eventWeight;
      if (isam==1) ttCount+=eventWeight;
      if (isam==2) bbCount+=eventWeight;
      if (isam==3) llCount+=eventWeight;
    
    } // end entry loop
    
    delete infile;
    infile=0, intree=0;

  } // end sample loop

  cout << "             HH \t tt \t BB \t LL " << endl;
  cout << "ntuple     : " << hhTotal << "\t" << ttTotal << "\t" << bbTotal << "\t" << llTotal << endl;
  cout << "preselect  : " << hhPre << "\t" << ttPre << "\t" << bbPre << "\t" << llPre << endl;
  cout << "m(tautau)  : " << hh1 << "\t" << tt1 << "\t" << bb1 << "\t" << ll1 << endl;
  cout << "m(bb)      : " << hh2 << "\t" << tt2 << "\t" << bb2 << "\t" << ll2 << endl;
  cout << "m(bbtautau): " << hh3 << "\t" << tt3 << "\t" << bb3 << "\t" << ll3 << endl;
  cout << "mT2        : " << hhCount << "\t" << ttCount << "\t" << bbCount << "\t" << llCount << endl;
  cout << "at 3/ab    : " << hhCount*3000 << "\t" << ttCount*3000 << "\t" << bbCount*3000 << "\t" << llCount*3000 << endl;
  cout << "S/sqrt(B)  : " << hhCount*3000/TMath::Sqrt(3000*(ttCount+bbCount+llCount)) << endl;

  Float_t L = 3000;
  Float_t hhxsec = 2.92;

  Float_t scale = 40.0/hhxsec;

  Float_t N = hhCount*L;
  Float_t dN = TMath::Sqrt(hhCount*L);

  Float_t Ae = hhCount/hhxsec;
  Float_t dAe = TMath::Sqrt(hhCount)/hhxsec;

  Float_t xsec = N/(Ae*L);
  Float_t dxsec = TMath::Sqrt(dN*dN/((Ae*L)*(Ae*L))+ dAe*dAe*N*N/(((Ae*L)*(Ae*L))));

  cout << xsec << " +- " << dxsec << " fb" << endl;
  cout << xsec*scale << " +- " << dxsec*scale << " fb" << endl;

  //p0    : 68.5935 +/- 0.187739
  //p1    : -44.6022 +/- 0.0552735
  //p2    : 9.13082 +/- 0.0204945

  TCanvas *c = MakeCanvas("c", "c", 600, 600);

  TF1 *smcurve = new TF1("smcurve", "[0]*(68.5935-44.6022*x+9.13082*x*x)",0,5);
  smcurve->SetParameter(0,1.0);
  cout << smcurve->Eval(1.0) << endl;
  smcurve->SetParameter(0,40/smcurve->Eval(1.0));
  cout << smcurve->Eval(1.0) << endl;

  smcurve->SetTitle("");
  smcurve->GetXaxis()->SetTitle("#lambda");
  smcurve->GetYaxis()->SetTitle("#sigma");
  
  smcurve->SetLineColor(kBlue);
  smcurve->Draw();
  
  TF1 *fxsec = new TF1("fxsec", "[0]", 0, 5);
  fxsec->SetParameter(0,xsec*scale+2*dxsec*scale);
  fxsec->SetLineColor(kRed);
  fxsec->Draw("same");

  cout << smcurve->GetX((xsec+2*dxsec)*scale,1,5) << endl;
  cout << smcurve->GetX((xsec+2*dxsec)*scale,0,1) << endl;

  /*    
  Float_t scale=1;

  CPlot plotTpt1("t_pt1", "", "Leading #tau p_{T} [GeV/c]","Events");
  for (Int_t i=0; i<sampleNames.size(); i++) {
    scale=hTpt1[i]->Integral();
    hTpt1[i]->Scale(1/scale);
    plotTpt1.AddHist1D(hTpt1[i], sampleTitles[i],"hist",sampleColors[i]);
  }
  plotTpt1.TransLegend(0.0,0);
  plotTpt1.Draw(c, kTRUE,"png",1);

  CPlot plotTpt2("t_pt2", "", "Second #tau p_{T} [GeV/c]","Events");
  for (Int_t i=0; i<sampleNames.size(); i++) {
    scale=hTpt2[i]->Integral();
    hTpt2[i]->Scale(1/scale);
    plotTpt2.AddHist1D(hTpt2[i], sampleTitles[i],"hist",sampleColors[i]);
  }
  plotTpt2.TransLegend(0.0,0);
  plotTpt2.Draw(c, kTRUE,"png",1);

  CPlot plotBpt1("b_pt1", "", "Leading b p_{T} [GeV/c]","Events");
  for (Int_t i=0; i<sampleNames.size(); i++) {
    scale=hBpt1[i]->Integral();
    hBpt1[i]->Scale(1/scale);
    plotBpt1.AddHist1D(hBpt1[i], sampleTitles[i],"hist",sampleColors[i]);
  }
  plotBpt1.TransLegend(0.0,0);
  plotBpt1.Draw(c, kTRUE,"png",1);

  CPlot plotBpt2("b_pt2", "", "Second b p_{T} [GeV/c]","Events");
  for (Int_t i=0; i<sampleNames.size(); i++) {
    scale=hBpt2[i]->Integral();
    hBpt2[i]->Scale(1/scale);
    plotBpt2.AddHist1D(hBpt2[i], sampleTitles[i],"hist",sampleColors[i]);
  }
  plotBpt2.TransLegend(0.0,0);
  plotBpt2.Draw(c, kTRUE,"png",1);

  CPlot plotTeta1("t_eta1", "", "Leading #tau |#eta|","Events");
  for (Int_t i=0; i<sampleNames.size(); i++) {
    scale=hTeta1[i]->Integral();
    hTeta1[i]->Scale(1/scale);
    plotTeta1.AddHist1D(hTeta1[i], sampleTitles[i],"hist",sampleColors[i]);
  }
  plotTeta1.TransLegend(0.0,0);
  plotTeta1.Draw(c, kTRUE,"png",1);

  CPlot plotTeta2("t_eta2", "", "Second #tau |#eta|","Events");
  for (Int_t i=0; i<sampleNames.size(); i++) {
    scale=hTeta2[i]->Integral();
    hTeta2[i]->Scale(1/scale);
    plotTeta2.AddHist1D(hTeta2[i], sampleTitles[i],"hist",sampleColors[i]);
  }
  plotTeta2.TransLegend(0.0,0);
  plotTeta2.Draw(c, kTRUE,"png",1);

  CPlot plotBeta1("b_eta1", "", "Leading b |#eta|","Events");
  for (Int_t i=0; i<sampleNames.size(); i++) {
    scale=hBeta1[i]->Integral();
    hBeta1[i]->Scale(1/scale);
    plotBeta1.AddHist1D(hBeta1[i], sampleTitles[i],"hist",sampleColors[i]);
  }
  plotBeta1.TransLegend(0.0,0);
  plotBeta1.Draw(c, kTRUE,"png",1);

  CPlot plotBeta2("b_eta2", "", "Second b |#eta|","Events");
  for (Int_t i=0; i<sampleNames.size(); i++) {
    scale=hBeta2[i]->Integral();
    hBeta2[i]->Scale(1/scale);
    plotBeta2.AddHist1D(hBeta2[i], sampleTitles[i],"hist",sampleColors[i]);
  }
  plotBeta2.TransLegend(0.0,0);
  plotBeta2.Draw(c, kTRUE,"png",1);

  CPlot plotHtt("m_tt", "", "M(#tau#tau) [GeV/c^{2}]","Events");
  for (Int_t i=0; i<sampleNames.size(); i++) {
    plotHtt.AddHist1D(hTTMass[i], sampleTitles[i],"hist",sampleColors[i]);
  }
  plotHtt.TransLegend(0.0,0);
  plotHtt.Draw(c, kTRUE,"png",1);

  CPlot plotHbb("m_bb", "", "M(bb) [GeV/c^{2}]","Events");
  for (Int_t i=0; i<sampleNames.size(); i++) {
    plotHbb.AddHist1D(hBBMass[i], sampleTitles[i],"hist",sampleColors[i]);
  }
  plotHbb.TransLegend(0.0,0);
  plotHbb.Draw(c, kTRUE, "png", 1);

  CPlot plotHttbb("m_ttbb", "", "M(#tau#taubb) [GeV/c^{2}]","Events");
  for (Int_t i=0; i<sampleNames.size(); i++) {
    plotHttbb.AddHist1D(hTTBBMass[i], sampleTitles[i],"hist",sampleColors[i]);
  }
  plotHttbb.TransLegend(0.0,0);
  plotHttbb.Draw(c, kTRUE,"png",1);

  CPlot plotMt2("mt2", "", "mT2","Events");
  for (Int_t i=0; i<sampleNames.size(); i++) {
    plotMt2.AddHist1D(hMt2[i], sampleTitles[i],"hist",sampleColors[i]);
  }
  plotMt2.TransLegend(0.0,0);
  plotMt2.Draw(c, kTRUE,"png",1);

  CPlot plotTTangle("ang_tt", "", "#Delta#phi(#tau#tau)","Events");
  for (Int_t i=0; i<sampleNames.size(); i++) {
    plotTTangle.AddHist1D(hTTAngle[i], sampleTitles[i],"hist",sampleColors[i]);
  }
  plotTTangle.TransLegend(0.0,0);
  plotTTangle.Draw(c, kTRUE,"png",1);
  
  CPlot plotBBangle("ang_bb", "", "#Delta#phi(bb)","Events");
  for (Int_t i=0; i<sampleNames.size(); i++) {
    plotBBangle.AddHist1D(hBBAngle[i], sampleTitles[i],"hist",sampleColors[i]);
  }
  plotBBangle.TransLegend(0.0,0);
  plotBBangle.Draw(c, kTRUE,"png",1);
  */
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
