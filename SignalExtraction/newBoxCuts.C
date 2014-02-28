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

void newBoxCuts(const TString conf="new.conf", Int_t tauType=0) {

  TGaxis::SetMaxDigits(3);

  // tau decay modes
  enum { hadron=1, electron, muon };

  // tautau decay modes
  enum { all=0, dijet, jetmu, jetele, muele };

  // background processes
  enum { HH=0, TT, ZH, WH, WW, ZZ, ZW, ETC };

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

  // read in file names
  vector<TString> sampleNames;
  vector<TString> sampleTitles;
  vector<Int_t> sampleColors;
  confParse(conf, sampleNames, sampleTitles, sampleColors);

  // define pre-selection variables
  const Float_t TAU_PT_MIN_H = 30;
  const Float_t TAU_PT_MIN_L = 20;
  const Float_t B_PT_MIN = 30;
  
  const Float_t TAU_ETA_MAX = 2.5;
  const Float_t B_ETA_MAX = 2.5;
  const Int_t ETA_BINS = 12;

  // signal/background processes
  vector<TString> processType;
  processType.push_back("HH"); processType.push_back("tt"); processType.push_back("ZH"); processType.push_back("WH"); 
  processType.push_back("WW"); processType.push_back("ZZ"); processType.push_back("ZW"); processType.push_back("ETC");

  vector<Double_t> totalXsec;
  vector<Double_t> preselectXsec;
  vector<Double_t> kineCutsXsec;
  vector<Double_t> tautauCutsXsec;
  vector<Double_t> bbCutsXsec;
  vector<Double_t> hhCutsXsec;
  vector<Double_t> finalXsec;

  vector<Double_t> preselectXsecUnc;
  vector<Double_t> kineCutsXsecUnc;
  vector<Double_t> bbCutsXsecUnc;
  vector<Double_t> hhCutsXsecUnc;
  vector<Double_t> finalXsecUnc;

  vector<TH1D*> hTpt1, hTpt2, hBpt1, hBpt2;
  vector<TH1D*> hTeta1, hTeta2, hBeta1, hBeta2, hBBdR;
  vector<TH1D*> hTTMass, hBBMass, hTTBBMass;
  vector<TH1D*> hTTAngle, hBBAngle;
  vector<TH1D*> hMt2;

  char hname[100];
  
  for(UInt_t isam=0; isam<processType.size(); isam++) {
    totalXsec.push_back(0); preselectXsec.push_back(0); kineCutsXsec.push_back(0); tautauCutsXsec.push_back(0); 
    bbCutsXsec.push_back(0); hhCutsXsec.push_back(0);    finalXsec.push_back(0);
    preselectXsecUnc.push_back(0); kineCutsXsecUnc.push_back(0);
    bbCutsXsecUnc.push_back(0); hhCutsXsecUnc.push_back(0);    finalXsecUnc.push_back(0);

    sprintf(hname, "hTpt1_%i",isam); hTpt1.push_back(new TH1D(hname, processType[isam], 100, 0, 300)); hTpt1[isam]->Sumw2(); 
    sprintf(hname, "hTpt2_%i",isam); hTpt2.push_back(new TH1D(hname, processType[isam], 100, 0, 300)); hTpt2[isam]->Sumw2(); 
    sprintf(hname, "hBpt1_%i",isam); hBpt1.push_back(new TH1D(hname, processType[isam], 100, 0, 300)); hBpt1[isam]->Sumw2(); 
    sprintf(hname, "hBpt2_%i",isam); hBpt2.push_back(new TH1D(hname, processType[isam], 100, 0, 300)); hBpt2[isam]->Sumw2(); 

    sprintf(hname, "hBBdR_%i",isam); hBBdR.push_back(new TH1D(hname, processType[isam], 16, 0, 8.0)); hBBdR[isam]->Sumw2(); 

    sprintf(hname, "hTeta1_%i",isam); hTeta1.push_back(new TH1D(hname, processType[isam], 16, -2.5, 2.5)); hTeta1[isam]->Sumw2(); 
    sprintf(hname, "hTeta2_%i",isam); hTeta2.push_back(new TH1D(hname, processType[isam], 16, -2.5, 2.5)); hTeta2[isam]->Sumw2(); 
    sprintf(hname, "hBeta1_%i",isam); hBeta1.push_back(new TH1D(hname, processType[isam], 16, -2.5, 2.5)); hBeta1[isam]->Sumw2(); 
    sprintf(hname, "hBeta2_%i",isam); hBeta2.push_back(new TH1D(hname, processType[isam], 16, -2.5, 2.5)); hBeta2[isam]->Sumw2(); 

    sprintf(hname, "hTTMass_%i",isam); hTTMass.push_back(new TH1D(hname, processType[isam], 100, 0, 400)); hTTMass[isam]->Sumw2(); 
    sprintf(hname, "hBBMass_%i",isam); hBBMass.push_back(new TH1D(hname, processType[isam], 100, 0, 400)); hBBMass[isam]->Sumw2(); 
    sprintf(hname, "hTTBBMass_%i",isam); hTTBBMass.push_back(new TH1D(hname, processType[isam], 100, 0, 1200)); hTTBBMass[isam]->Sumw2(); 
    sprintf(hname, "hTTAngle_%i",isam); hTTAngle.push_back(new TH1D(hname, processType[isam], 20, 0, TMath::Pi())); hTTAngle[isam]->Sumw2(); 
    sprintf(hname, "hBBAngle_%i",isam); hBBAngle.push_back(new TH1D(hname, processType[isam], 20, 0, TMath::Pi())); hBBAngle[isam]->Sumw2(); 
    sprintf(hname, "hMt2_%i",isam); hMt2.push_back(new TH1D(hname, processType[isam], 50, 0, 500)); hMt2[isam]->Sumw2(); 
  }

  Float_t eventWeight=1;
  UInt_t eventType;
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
    intree->SetBranchAddress("recoB1",         &recoB1);       // 4-vector for reconstructed leading b-jet
    intree->SetBranchAddress("recoB2",         &recoB2);       // 4-vector for reconstructed second b-jet
    intree->SetBranchAddress("recoExtraJet",   &recoExtraJet); // 4-vector for reconstructed extra jet

    for(UInt_t iEntry=0; iEntry<intree->GetEntries(); iEntry++) { // entry loop
      intree->GetEntry(iEntry);

      totalXsec[eventType]+=eventWeight;

      // skip events that don't have 2 reco b's and 2 reco tau's
      if ( (recoB1->Pt()==999) || (recoB2->Pt()==999) ) continue;
      if ( (recoTau1->Pt()==999) || (recoTau2->Pt()==999) ) continue;

      // exclude ditau->ee/mumu channels
      if ( (tauCat1==electron) && (tauCat2==electron) ) continue;
      if ( (tauCat1==muon) && (tauCat2==muon) ) continue;

      // select tautypes

      if ( (tauType==dijet) && ( ( tauCat1 != hadron ) || ( tauCat2!=hadron ) ) ) continue;

      if ( (tauType==jetmu) && !( ( ( tauCat1==hadron ) && ( tauCat2==muon) ) || ( ( tauCat2==hadron ) && ( tauCat1==muon) ) ) ) continue;

      if ( (tauType==jetele) && !( ( (tauCat1==hadron) && (tauCat2==electron) ) || ( (tauCat2==hadron) && (tauCat1==electron) ) ) ) continue;

      if ( (tauType==muele) && !( ( (tauCat1==muon) && (tauCat2==electron) ) || ( (tauCat2==muon) && (tauCat1==electron) ) ) ) continue;

      bPt1=recoB1->Pt(); 
      bPt2=recoB2->Pt();

      if (tauType==jetmu) {
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
      if (tauType==jetele) {
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
      if (tauType==muele) {
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

      preselectXsec[eventType]+=eventWeight;
      preselectXsecUnc[eventType]+=eventWeight*eventWeight;

      hBBdR[eventType]->Fill(deltaR(recoB1->Eta(), recoB2->Eta(), recoB1->Phi(), recoB2->Phi()), eventWeight);

      // tau kinematics
      if ( ( fabs( recoTau1->Eta() ) > TAU_ETA_MAX ) || ( fabs(recoTau2->Eta() ) > TAU_ETA_MAX ) ) continue;

      if (( tauCat1 == hadron ) && (tauPt1 < TAU_PT_MIN_H )) continue;
      else if (( tauCat1 != hadron ) && (tauPt1 < TAU_PT_MIN_L )) continue;

      if (( tauCat2 == hadron ) && (tauPt2 < TAU_PT_MIN_H )) continue;
      else if (( tauCat2 != hadron ) && (tauPt2 < TAU_PT_MIN_L )) continue;

      // b kinematics
      if ( ( fabs( recoB1->Eta() ) > B_ETA_MAX ) || ( fabs( recoB2->Eta() ) > B_ETA_MAX ) ) continue;

      if ( ( bPt1 < B_PT_MIN ) || ( bPt2 < B_PT_MIN ) ) continue;

      kineCutsXsec[eventType]+=eventWeight;
      kineCutsXsecUnc[eventType]+=eventWeight*eventWeight;

      // require boost of at least 50 GeV
      //if ( (recoExtraJet->Pt()!=999) && (recoExtraJet->Pt()>50) ) continue;

      hTpt1[eventType]->Fill(tauPt1, eventWeight);
      hTpt2[eventType]->Fill(tauPt2, eventWeight);
      hBpt1[eventType]->Fill(bPt1, eventWeight);
      hBpt2[eventType]->Fill(bPt2, eventWeight);

      hTeta1[eventType]->Fill(recoTau1->Eta(), eventWeight);
      hTeta2[eventType]->Fill(recoTau2->Eta(), eventWeight);
      hBeta1[eventType]->Fill(recoB1->Eta(), eventWeight);
      hBeta2[eventType]->Fill(recoB2->Eta(), eventWeight);

      LorentzVector vTau1(tauPt1, recoTau1->Eta(), recoTau1->Phi(), recoTau1->M());
      LorentzVector vTau2(tauPt2, recoTau2->Eta(), recoTau2->Phi(), recoTau2->M());
      LorentzVector vTauHiggs = vTau1 + vTau2;
      
      LorentzVector vB1(bPt1, recoB1->Eta(), recoB1->Phi(), recoB1->M());
      LorentzVector vB2(bPt2, recoB2->Eta(), recoB2->Phi(), recoB2->M());
      LorentzVector vBHiggs = vB1 + vB2;

      LorentzVector vHH = vTauHiggs + vBHiggs;

      for (Int_t c=0; c<6; c++) {

	if ( (c!=0) && ((vTauHiggs.M()<50) || (vTauHiggs.M()>120)) ) continue;
	if ( (c!=1) && ((vBHiggs.M()<100) || (vBHiggs.M()>150)) ) continue;
	if ( (c!=2) && (vHH.M()<360) ) continue;
	if ( (c!=3) && (mt2<120) ) continue;

	if (c==0) hTTMass[eventType]->Fill(vTauHiggs.M(),eventWeight);
	else if (c==1) hBBMass[eventType]->Fill(vBHiggs.M(),eventWeight);
	else if (c==2) hTTBBMass[eventType]->Fill(vHH.M(),eventWeight);
	else if (c==3) hMt2[eventType]->Fill(mt2,eventWeight);
	else if (c==4) hTTAngle[eventType]->Fill(deltaPhi(recoTau1->Phi(), recoTau2->Phi()),eventWeight);
	else if (c==5) hBBAngle[eventType]->Fill(deltaPhi(recoB1->Phi(), recoB2->Phi()),eventWeight);

      }

      if ( (vTauHiggs.M()<50) || (vTauHiggs.M()>120) ) continue;

      tautauCutsXsec[eventType]+=eventWeight;

      if ( (vBHiggs.M()<100) || (vBHiggs.M()>150) ) continue;

      bbCutsXsec[eventType]+=eventWeight;
      bbCutsXsecUnc[eventType]+=eventWeight*eventWeight;
      
      if ( vHH.M()<360 ) continue;

      hhCutsXsec[eventType]+=eventWeight;
      hhCutsXsecUnc[eventType]+=eventWeight*eventWeight;
      
      if ( mt2<120 ) continue;
      
      finalXsec[eventType]+=eventWeight;
      finalXsecUnc[eventType]+=eventWeight*eventWeight;
    } // end entry loop
    
    delete infile;
    infile=0, intree=0;

  } // end sample loop

  vector<Double_t> totalBackgrounds;
  for (Int_t j=0; j<7; j++) {
    totalBackgrounds.push_back(0);
  }

  cout << "sample\t\t presel\t tau/b\t bb\t hh\t final" << endl;

  for (UInt_t i=0; i<processType.size(); i++) {
    cout << processType[i] << "\t" << preselectXsec[i]*3000 << "\t" <<  kineCutsXsec[i]*3000 << "\t" << bbCutsXsec[i]*3000 << "\t";
    cout << hhCutsXsec[i]*3000 << "\t" << finalXsec[i]*3000 << endl;
    cout << processType[i] << "\t" << TMath::Sqrt(preselectXsecUnc[i])*3000 << "\t" <<  TMath::Sqrt(kineCutsXsecUnc[i])*3000 << "\t" << TMath::Sqrt(bbCutsXsecUnc[i])*3000 << "\t";
    cout << TMath::Sqrt(hhCutsXsecUnc[i])*3000 << "\t" << TMath::Sqrt(finalXsecUnc[i])*3000 << endl;

    if (i!=0) {
      totalBackgrounds[0]+=totalXsec[i]; totalBackgrounds[1]+=preselectXsec[i]; totalBackgrounds[2]+=kineCutsXsec[i]; totalBackgrounds[3]+=tautauCutsXsec[i]; 
      totalBackgrounds[4]+=bbCutsXsec[i]; totalBackgrounds[5]+=hhCutsXsec[i]; totalBackgrounds[6]+=finalXsec[i];
    }
  }

  cout << "S/sqrt(B): " << TMath::Sqrt(3000)*totalXsec[0]/TMath::Sqrt(totalBackgrounds[0]) << "\t" << TMath::Sqrt(3000)*kineCutsXsec[0]/TMath::Sqrt(totalBackgrounds[2]) << "\t";
  cout << TMath::Sqrt(3000)*bbCutsXsec[0]/TMath::Sqrt(totalBackgrounds[4]) << "\t" << TMath::Sqrt(3000)*hhCutsXsec[0]/TMath::Sqrt(totalBackgrounds[5]) << "\t" << TMath::Sqrt(3000)*finalXsec[0]/TMath::Sqrt(totalBackgrounds[6]) << endl;

  /*
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
  */
  TCanvas *c = MakeCanvas("c", "c", 600, 600);
  /*
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
  */

  Float_t scale=1;

  CPlot plotBBdR("bb_dr", "", "#Deltar(bb)","Events");
  for (Int_t i=0; i<sampleNames.size(); i++) {
    scale=hBBdR[i]->Integral();
    hBBdR[i]->Scale(1/scale);
    plotBBdR.AddHist1D(hBBdR[i], sampleTitles[i],"hist",sampleColors[i]);
  }
  plotBBdR.TransLegend(0.0,0);
  plotBBdR.Draw(c, kTRUE,"png",1);

  /*
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
    scale=hTTMass[i]->Integral();
    hTTMass[i]->Scale(1/scale);
    plotHtt.AddHist1D(hTTMass[i], sampleTitles[i],"hist",sampleColors[i]);
  }
  plotHtt.TransLegend(0.0,0);
  plotHtt.Draw(c, kTRUE,"png",1);

  CPlot plotHbb("m_bb", "", "M(bb) [GeV/c^{2}]","Events");
  for (Int_t i=0; i<sampleNames.size(); i++) {
    scale=hBBMass[i]->Integral();
    hBBMass[i]->Scale(1/scale);
    plotHbb.AddHist1D(hBBMass[i], sampleTitles[i],"hist",sampleColors[i]);
  }
  plotHbb.TransLegend(0.0,0);
  plotHbb.Draw(c, kTRUE, "png", 1);

  CPlot plotHttbb("m_ttbb", "", "M(#tau#taubb) [GeV/c^{2}]","Events");
  for (Int_t i=0; i<sampleNames.size(); i++) {
    scale=hTTBBMass[i]->Integral();
    hTTBBMass[i]->Scale(1/scale);
    plotHttbb.AddHist1D(hTTBBMass[i], sampleTitles[i],"hist",sampleColors[i]);
  }
  plotHttbb.TransLegend(0.0,0);
  plotHttbb.Draw(c, kTRUE,"png",1);

  CPlot plotMt2("mt2", "", "mT2","Events");
  for (Int_t i=0; i<sampleNames.size(); i++) {
    scale=hMt2[i]->Integral();
    hMt2[i]->Scale(1/scale);
    plotMt2.AddHist1D(hMt2[i], sampleTitles[i],"hist",sampleColors[i]);
  }
  plotMt2.TransLegend(0.0,0);
  plotMt2.Draw(c, kTRUE,"png",1);

  CPlot plotTTangle("ang_tt", "", "#Delta#phi(#tau#tau)","Events");
  for (Int_t i=0; i<sampleNames.size(); i++) {
    scale=hTTAngle[i]->Integral();
    hTTAngle[i]->Scale(1/scale);
    plotTTangle.AddHist1D(hTTAngle[i], sampleTitles[i],"hist",sampleColors[i]);
  }
  plotTTangle.TransLegend(0.0,0);
  plotTTangle.Draw(c, kTRUE,"png",1);
  
  CPlot plotBBangle("ang_bb", "", "#Delta#phi(bb)","Events");
  for (Int_t i=0; i<sampleNames.size(); i++) {
    scale=hBBAngle[i]->Integral();
    hBBAngle[i]->Scale(1/scale);
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
