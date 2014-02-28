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
#include <TLorentzVector.h>
#include "Math/LorentzVector.h"

#include "modules/Delphes.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "classes/DelphesClasses.h"

#include "../MitStyleRemix.hh"
#include "../CPlot.hh"

#endif

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;

void confParse(const TString conf, vector<TString> &sampleNames, vector<TString> &sampleTitles, vector<Int_t> &sampleColors);

Float_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 );
Double_t deltaPhi(const Double_t phi1, const Double_t phi2);

void distributions(const TString conf="new.conf") {

  CPlot::sOutDir = "noboost50gev"; 

  // define kinematic/plotting constants
  const Float_t TAU_PT_MIN_H = 30;
  const Float_t TAU_PT_MIN_L = 20;
  const Float_t B_PT_MIN = 30;
  
  const Float_t TAU_ETA_MAX = 2.5;
  const Float_t B_ETA_MAX = 2.5;
  const Int_t ETA_BINS = 12;

  const Float_t TAU_PT_MAX = 500;
  const Float_t B_PT_MAX = 500;

  const Float_t HTT_MAX = 180;
  const Float_t HTT_MIN = 60;
  const Float_t HBB_MAX = 160;
  const Float_t HBB_MIN = 80;

  // tau decay modes
  enum { hadron=1, electron, muon };
  // tau-tau decay modes
  enum {dijet=0, jetele, jetmu, elemu };

  vector<TString> sampleNames;
  vector<TString> sampleTitles;
  vector<Int_t> sampleColors;

  confParse(conf, sampleNames, sampleTitles, sampleColors);

  //cat: jet-jet, jet-ele, jet-mu, ele-mu
  vector<TH1D*> vTTMass[4], vBBMass[4], vTTBBMass[4];
  vector<TH1D*> vTTAngle[4], vBBAngle[4], vTTBBAngle[4];
  vector<TH1D*> vHPt[4];
  vector<TH1D*> hTTMass, hBBMass, hTTBBMass;
  vector<TH1D*> hTTAngle, hBBAngle, hTTBBAngle;
  vector<TH1D*> hHPt;

  char hname[100];
  
  for (UInt_t icat=0; icat<4; icat++) {
    for(UInt_t isam=0; isam<sampleNames.size(); isam++) {
      sprintf(hname, "hTTMass_cat%i_%i",icat,isam);   vTTMass[icat].push_back(new TH1D(hname, "", 50, 0, 600));     vTTMass[icat][isam]->Sumw2();
      sprintf(hname, "hBBMass_cat%i_%i",icat,isam);   vBBMass[icat].push_back(new TH1D(hname, "", 50, 0, 600));     vBBMass[icat][isam]->Sumw2();
      sprintf(hname, "hTTBBMass_cat%i_%i",icat,isam); vTTBBMass[icat].push_back(new TH1D(hname, "", 50, 0, 1200)); vTTBBMass[icat][isam]->Sumw2();

      sprintf(hname, "hTTAngle_cat%i_%i",icat,isam);   vTTAngle[icat].push_back(new TH1D(hname, "", 20, 0, TMath::Pi()));     vTTAngle[icat][isam]->Sumw2();
      sprintf(hname, "hBBAngle_cat%i_%i",icat,isam);   vBBAngle[icat].push_back(new TH1D(hname, "", 20, 0, TMath::Pi()));     vBBAngle[icat][isam]->Sumw2();
      sprintf(hname, "hTTBBAngle_cat%i_%i",icat,isam); vTTBBAngle[icat].push_back(new TH1D(hname, "", 20, 0, TMath::Pi()));   vTTBBAngle[icat][isam]->Sumw2();

      sprintf(hname, "hHPt_cat%i_%i",icat,isam); vHPt[icat].push_back(new TH1D(hname, "", 50, 0, 500)); vHPt[icat][isam]->Sumw2();

    }
  }

  for(UInt_t isam=0; isam<sampleNames.size(); isam++) {
    sprintf(hname, "hTTMass_%i",isam); hTTMass.push_back(new TH1D(hname, sampleTitles[isam], 100, 0, 600)); hTTMass[isam]->Sumw2(); 
    sprintf(hname, "hBBMass_%i",isam); hBBMass.push_back(new TH1D(hname, sampleTitles[isam], 100, 0, 600)); hBBMass[isam]->Sumw2(); 
    sprintf(hname, "hTTBBMass_%i",isam); hTTBBMass.push_back(new TH1D(hname, sampleTitles[isam], 100, 0, 1200)); hTTBBMass[isam]->Sumw2(); 
    sprintf(hname, "hTTAngle_%i",isam); hTTAngle.push_back(new TH1D(hname, sampleTitles[isam], 20, 0, TMath::Pi())); hTTAngle[isam]->Sumw2(); 
    sprintf(hname, "hBBAngle_%i",isam); hBBAngle.push_back(new TH1D(hname, sampleTitles[isam], 20, 0, TMath::Pi())); hBBAngle[isam]->Sumw2(); 
    sprintf(hname, "hTTBBAngle_%i",isam); hTTBBAngle.push_back(new TH1D(hname, sampleTitles[isam], 20, 0, TMath::Pi())); hTTBBAngle[isam]->Sumw2(); 
    sprintf(hname, "hHPt_%i",isam); hHPt.push_back(new TH1D(hname, sampleTitles[isam], 100, 0, 500)); hHPt[isam]->Sumw2(); 
  }

  //cout << endl;
  //cout << " --- Applied Cuts: --- " << endl;
  //cout << "b   pT > " << B_PT_MIN << " and |eta| < " << B_ETA_MAX << endl;
  //cout << "tau pT > " << TAU_PT_MIN << " and |eta| < " << TAU_ETA_MAX << endl;
  //cout << endl;

  Float_t eventWeight=1;
  UInt_t bTag1, bTag2;
  UInt_t tauDecayCat1, tauDecayCat2;
  LorentzVector *recoB1=0, *recoB2=0;
  LorentzVector *recoTau1=0, *recoTau2=0;
  LorentzVector *recoExtraJet=0;
  LorentzVector *tauHiggs=0, *bHiggs=0;

  TFile *infile;
  TTree *intree;

  Float_t bPt1=0, bPt2=0;
  Float_t tauPt1=0, tauPt2=0;

  for (UInt_t isam=0; isam<sampleNames.size(); isam++) { // sample loop

    TString infilename = sampleNames[isam];
    //cout << "Processing  " << infilename << " ..." << endl;
    infile = new TFile(infilename); assert(infile);
    intree = (TTree*) infile->Get("Events"); assert(intree);
 
    intree->SetBranchAddress("eventWeight",    &eventWeight);
    intree->SetBranchAddress("bTag1",          &bTag1);
    intree->SetBranchAddress("bTag2",          &bTag2);
    //intree->SetBranchAddress("genB1",          &genB1);
    //intree->SetBranchAddress("genB2",          &genB2);
    intree->SetBranchAddress("recoB1",         &recoB1);
    intree->SetBranchAddress("recoB2",         &recoB2);
    intree->SetBranchAddress("tauCat1",   &tauDecayCat1);
    intree->SetBranchAddress("tauCat2",   &tauDecayCat2);
    //intree->SetBranchAddress("genTau1",        &genTau1);
    //intree->SetBranchAddress("genTau2",        &genTau2);
    //intree->SetBranchAddress("genDecayTau1",   &genDecayTau1);
    // intree->SetBranchAddress("genDecayTau2",   &genDecayTau2);
    intree->SetBranchAddress("recoTau1",       &recoTau1);
    intree->SetBranchAddress("recoTau2",       &recoTau2);
    intree->SetBranchAddress("recoExtraJet",   &recoExtraJet); 

    for(UInt_t iEntry=0; iEntry<intree->GetEntries(); iEntry++) { // entry loop
      intree->GetEntry(iEntry);

      // skip events that don't have 2 reco b's and 2 reco tau's
      if ( (recoB1->Pt()==999) || (recoB2->Pt()==999) ) continue;
      if ( (recoTau1->Pt()==999) || (recoTau2->Pt()==999) ) continue;

      if ( (recoExtraJet->Pt()!=999) && (recoExtraJet->Pt()>50) ) continue;

      bPt1=recoB1->Pt(); 
      bPt2=recoB2->Pt();
      tauPt1 = recoTau1->Pt();
      tauPt2 = recoTau2->Pt();

      if ((tauDecayCat1 == electron) && (tauDecayCat2 == electron)) continue;
      if ((tauDecayCat1 == muon) && (tauDecayCat2 == muon)) continue;

      // tau pre-selection
      if ( ( fabs( recoTau1->Eta() ) > TAU_ETA_MAX ) || ( fabs(recoTau2->Eta() ) > TAU_ETA_MAX ) ) continue;
      if (( tauDecayCat1 == hadron ) && (tauPt1 < TAU_PT_MIN_H )) continue;
      else if (( tauDecayCat2 == hadron ) && (tauPt2 < TAU_PT_MIN_H )) continue;
      else if ( ( tauPt1 < TAU_PT_MIN_L ) || ( tauPt2 < TAU_PT_MIN_L ) ) continue;

      // b pre-selection
      if ( ( fabs( recoB1->Eta() ) > B_ETA_MAX ) || ( fabs( recoB2->Eta() ) > B_ETA_MAX ) ) continue;
      if ( ( bPt1 < B_PT_MIN ) || ( bPt2 < B_PT_MIN ) ) continue;

      LorentzVector vTau1(tauPt1, recoTau1->Eta(), recoTau1->Phi(), recoTau1->M());
      LorentzVector vTau2(tauPt2, recoTau2->Eta(), recoTau2->Phi(), recoTau2->M());
      LorentzVector vTauHiggs = vTau1 + vTau2;
      
      LorentzVector vB1(bPt1, recoB1->Eta(), recoB1->Phi(), recoB1->M());
      LorentzVector vB2(bPt2, recoB2->Eta(), recoB2->Phi(), recoB2->M());
      LorentzVector vBHiggs = vB1 + vB2;

      LorentzVector vHH = vTauHiggs + vBHiggs;

      hTTMass[isam]->Fill(vTauHiggs.M(),eventWeight);
      hBBMass[isam]->Fill(vBHiggs.M(),eventWeight);
      hTTBBMass[isam]->Fill(vHH.M(),eventWeight);

      hHPt[isam]->Fill(vHH.Pt(),eventWeight);

      hTTAngle[isam]->Fill(deltaPhi(recoTau1->Phi(), recoTau2->Phi()),eventWeight);
      hBBAngle[isam]->Fill(deltaPhi(recoB1->Phi(), recoB2->Phi()),eventWeight);
      hTTBBAngle[isam]->Fill(deltaPhi(vBHiggs.Phi(), vTauHiggs.Phi()),eventWeight);

      Int_t icat=-1;
      if ((tauDecayCat1 == hadron) && (tauDecayCat2 == hadron)) icat=dijet;
      else if ((tauDecayCat1 == hadron) && (tauDecayCat2 == electron)) icat=jetele;
      else if ((tauDecayCat1 == hadron) && (tauDecayCat2 == muon)) icat=jetmu;
      else if ((tauDecayCat1 == electron) && (tauDecayCat2 == muon)) icat=elemu;
      else if ((tauDecayCat1 == electron) && (tauDecayCat2 == hadron)) icat=jetele;
      else if ((tauDecayCat1 == muon) && (tauDecayCat2 == hadron)) icat=jetmu;
      else if ((tauDecayCat1 == muon) && (tauDecayCat2 == electron)) icat=elemu;
      if (icat<0) cout << "categories broke" << endl;

      vTTMass[icat][isam]->Fill(vTauHiggs.M(),eventWeight);
      vBBMass[icat][isam]->Fill(vBHiggs.M(),eventWeight);
      vTTBBMass[icat][isam]->Fill(vHH.M(),eventWeight);

      vHPt[icat][isam]->Fill(vHH.Pt(),eventWeight);

      vTTAngle[icat][isam]->Fill(deltaPhi(recoTau1->Phi(), recoTau2->Phi()),eventWeight);
      vBBAngle[icat][isam]->Fill(deltaPhi(recoB1->Phi(), recoB2->Phi()),eventWeight);
      vTTBBAngle[icat][isam]->Fill(deltaPhi(vBHiggs.Phi(), vTauHiggs.Phi()),eventWeight);

    } // end entry loop

    delete infile;
    infile=0, intree=0;

  } // end sample loop

  TCanvas *c = MakeCanvas("c", "c", 600, 600);

  Float_t scale=1;

  CPlot plotHtt("m_tt", "", "M(#tau#tau) [GeV/c^{2}]","Normalized Events");
  for (Int_t i=0; i<sampleNames.size(); i++) {
    scale=hTTMass[i]->Integral();
    hTTMass[i]->Scale(1/scale);    
    plotHtt.AddHist1D(hTTMass[i], sampleTitles[i],"hist",sampleColors[i]);
  }
  plotHtt.AddTextBox("All #tau channels",0.63,0.92,0.95,0.99,0);
  plotHtt.TransLegend(0.0,0);
  plotHtt.Draw(c, kTRUE,"png",1);

  CPlot plotHbb("m_bb", "", "M(bb) [GeV/c^{2}]","Normalized Events");
  for (Int_t i=0; i<sampleNames.size(); i++) {
    scale=hBBMass[i]->Integral();
    hBBMass[i]->Scale(1/scale);    
    plotHbb.AddHist1D(hBBMass[i], sampleTitles[i],"hist",sampleColors[i]);
  }
  plotHbb.AddTextBox("All #tau channels",0.63,0.92,0.95,0.99,0);
  plotHbb.TransLegend(0.0,0);
  plotHbb.Draw(c, kTRUE, "png", 1);

  CPlot plotHttbb("m_ttbb", "", "M(#tau#taubb) [GeV/c^{2}]","Normalized Events");
  for (Int_t i=0; i<sampleNames.size(); i++) {
    scale=hTTBBMass[i]->Integral();
    hTTBBMass[i]->Scale(1/scale);    
    plotHttbb.AddHist1D(hTTBBMass[i], sampleTitles[i],"hist",sampleColors[i]);
  }
  plotHttbb.AddTextBox("All #tau channels",0.63,0.92,0.95,0.99,0);
  plotHttbb.TransLegend(0.0,0);
  plotHttbb.Draw(c, kTRUE,"png",1);

  CPlot plotHPt("pt_ttbb", "", "p_{T}(#tau#taubb) [GeV/c]","Normalized Events");
  for (Int_t i=0; i<sampleNames.size(); i++) {
    scale=hHPt[i]->Integral();
    hHPt[i]->Scale(1/scale);    
    plotHPt.AddHist1D(hHPt[i], sampleTitles[i],"hist",sampleColors[i]);
  }
  plotHPt.AddTextBox("All #tau channels",0.63,0.92,0.95,0.99,0);
  plotHPt.TransLegend(0.0,0);
  plotHPt.Draw(c, kTRUE,"png",1);

  CPlot plotTTangle("ang_tt", "", "#Delta#phi(#tau#tau)","Normalized Events");
  for (Int_t i=0; i<sampleNames.size(); i++) {
    scale=hTTAngle[i]->Integral();
    hTTAngle[i]->Scale(1/scale);    
    plotTTangle.AddHist1D(hTTAngle[i], sampleTitles[i],"hist",sampleColors[i]);
  }
  plotTTangle.AddTextBox("All #tau channels",0.63,0.92,0.95,0.99,0);
  plotTTangle.TransLegend(0.0,0);
  plotTTangle.Draw(c, kTRUE,"png",1);
  
  CPlot plotBBangle("ang_bb", "", "#Delta#phi(bb)","Normalized Events");
  for (Int_t i=0; i<sampleNames.size(); i++) {
    scale=hBBAngle[i]->Integral();
    hBBAngle[i]->Scale(1/scale);    
    plotBBangle.AddHist1D(hBBAngle[i], sampleTitles[i],"hist",sampleColors[i]);
  }
  plotBBangle.AddTextBox("All #tau channels",0.63,0.92,0.95,0.99,0);
  plotBBangle.TransLegend(0.0,0);
  plotBBangle.Draw(c, kTRUE,"png",1);

  CPlot plotTTBBangle("ang_ttbb", "", "#Delta#phi(#tau#taubb)","Normalized Events");
  for (Int_t i=0; i<sampleNames.size(); i++) {
    scale=hTTBBAngle[i]->Integral();
    hTTBBAngle[i]->Scale(1/scale);    
    plotTTBBangle.AddHist1D(hTTBBAngle[i], sampleTitles[i],"hist",sampleColors[i]);
  }
  plotTTBBangle.AddTextBox("All #tau channels",0.63,0.92,0.95,0.99,0);
  plotTTBBangle.TransLegend(0.0,0);
  plotTTBBangle.Draw(c, kTRUE,"png",1);
  
  char title[100];

  for (Int_t i=0; i<4; i++) {

    if (i==dijet) sprintf(title,"Dijet #tau only");
    else if (i==jetele) sprintf(title,"Jet-ele #tau only");
    else if (i==jetmu) sprintf(title,"Jet-mu #tau only");
    else if (i==elemu) sprintf(title,"Ele-mu #tau only");

    sprintf(hname,"m_tt_%i",i);
    CPlot plotHtt2(hname, "", "M(#tau#tau) [GeV/c^{2}]","Normalized Events");
    for (Int_t j=0; j<sampleNames.size(); j++) {
      scale=vTTMass[i][j]->Integral();
      cout << title << ": " << sampleTitles[j] << ": " << scale << endl;
      vTTMass[i][j]->Scale(1/scale);    
      plotHtt2.AddHist1D(vTTMass[i][j], sampleTitles[j],"hist",sampleColors[j]);
    }
    plotHtt2.AddTextBox(title,0.63,0.92,0.95,0.99,0);
    plotHtt2.TransLegend(0.0,0);
    plotHtt2.Draw(c, kTRUE,"png",1);

    sprintf(hname,"m_bb_%i",i);
    CPlot plotHbb2(hname, "", "M(bb) [GeV/c^{2}]","Normalized Events");
    for (Int_t j=0; j<sampleNames.size(); j++) {
      scale=vBBMass[i][j]->Integral();
      vBBMass[i][j]->Scale(1/scale);    
      plotHbb2.AddHist1D(vBBMass[i][j], sampleTitles[j],"hist",sampleColors[j]);
    }
    plotHbb2.AddTextBox(title,0.63,0.92,0.95,0.99,0);
    plotHbb2.TransLegend(0.0,0);
    plotHbb2.Draw(c, kTRUE,"png",1);

    sprintf(hname,"m_ttbb_%i",i);
    CPlot plotHttbb2(hname, "", "M(#tau#taubb) [GeV/c^{2}]","Normalized Events");
    for (Int_t j=0; j<sampleNames.size(); j++) {
      scale=vTTBBMass[i][j]->Integral();
      vTTBBMass[i][j]->Scale(1/scale);    
      plotHttbb2.AddHist1D(vTTBBMass[i][j], sampleTitles[j],"hist",sampleColors[j]);
    }
    plotHttbb2.AddTextBox(title,0.63,0.92,0.95,0.99,0);
    plotHttbb2.TransLegend(0.0,0);
    plotHttbb2.Draw(c, kTRUE,"png",1);

    sprintf(hname,"pt_ttbb_%i",i);
    CPlot plotHPt2(hname, "", "p_{T}(#tau#taubb) [GeV/c]","Normalized Events");
    for (Int_t j=0; j<sampleNames.size(); j++) {
      scale=vHPt[i][j]->Integral();
      vHPt[i][j]->Scale(1/scale);    
      plotHPt2.AddHist1D(vHPt[i][j], sampleTitles[j],"hist",sampleColors[j]);
    }
    plotHPt2.AddTextBox(title,0.63,0.92,0.95,0.99,0);
    plotHPt2.TransLegend(0.0,0);
    plotHPt2.Draw(c, kTRUE,"png",1);

    sprintf(hname,"ang_tt_%i",i);
    CPlot plotTTangle2(hname, "", "#Delta#phi(#tau#tau)","Normalized Events");
    for (Int_t j=0; j<sampleNames.size(); j++) {
      scale=vTTAngle[i][j]->Integral();
      vTTAngle[i][j]->Scale(1/scale);    
      plotTTangle2.AddHist1D(vTTAngle[i][j], sampleTitles[j],"hist",sampleColors[j]);
    }
    plotTTangle2.AddTextBox(title,0.63,0.92,0.95,0.99,0);
    plotTTangle2.TransLegend(0.0,0);
    plotTTangle2.Draw(c, kTRUE,"png",1);
    
    sprintf(hname,"ang_bb_%i",i);
    CPlot plotBBangle2(hname, "", "#Delta#phi(bb)","Normalized Events");
    for (Int_t j=0; j<sampleNames.size(); j++) {
      scale=vBBAngle[i][j]->Integral();
      vBBAngle[i][j]->Scale(1/scale);    
      plotBBangle2.AddHist1D(vBBAngle[i][j], sampleTitles[j],"hist",sampleColors[j]);
    }
    plotBBangle2.AddTextBox(title,0.63,0.92,0.95,0.99,0);
    plotBBangle2.TransLegend(0.0,0);
    plotBBangle2.Draw(c, kTRUE,"png",1);

    sprintf(hname,"ang_ttbb_%i",i);
    CPlot plotTTBBangle2(hname, "", "#Delta#phi(#tau#taubb)","Normalized Events");
    for (Int_t j=0; j<sampleNames.size(); j++) {
      scale=vTTBBAngle[i][j]->Integral();
      vTTBBAngle[i][j]->Scale(1/scale);    
      plotTTBBangle2.AddHist1D(vTTBBAngle[i][j], sampleTitles[j],"hist",sampleColors[j]);
    }
    plotTTBBangle2.AddTextBox(title,0.63,0.92,0.95,0.99,0);
    plotTTBBangle2.TransLegend(0.0,0);
    plotTTBBangle2.Draw(c, kTRUE,"png",1);
    
  }

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
