#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include <TSystem.h>
#include <TNtuple.h>
#include <TH1D.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <TLorentzVector.h>
#include <assert.h>
#include <iomanip>
#endif

double DeltaPhi(double phi1, double phi2);

// root -l -q -b makeNtupleLHE_WZ.C+'("/afs/cern.ch/work/c/ceballos/public/samples/wzlhe8tev_qcdewk/","wzgamma_qed_5_qcd_99_sm.lhe",1,1,1)'
// root -l -q -b makeNtupleLHE_WZ.C+'("/afs/cern.ch/work/c/ceballos/public/samples/wwsslhe8tev_qcdewk/","qed_4_qcd_99_sm.lhe",1,1,0)'

//void makeNtupleLHE(TString pathDir="root://eoscms.cern.ch//store/group/phys_higgs/future/sixie/Madgraph/BBJJ_M60To200_14TeV_v2/LHE/",
void makeNtupleLHE(TString pathDir=".",
		   TString infname="events.99.lhe",
		   double weightIni = 1.0) {
  ifstream ifs(Form("%s/%s",pathDir.Data(),infname.Data()));
  assert(ifs.is_open());
  TString outNtuplename = Form("%s",infname.Data());
  outNtuplename.ReplaceAll(".lhe",".root");
  TFile *outtuple = TFile::Open(outNtuplename.Data(),"recreate");
  //TNtuple *nt = new TNtuple("Events","Events","ptq1:ptq2:ptq3:ptq4:etaq1:etaq2:etaq3:etaq4:phiq1:phiq2:phiq3:phiq4:pid1:pid2:pid3:pid4:weight");
  TNtuple *nt = new TNtuple("Events","Events","ptq1:ptq2:ptq3:ptq4:etaq1:etaq2:etaq3:etaq4:phiq1:phiq2:phiq3:phiq4:pid:weight");
  int eventType[10] = {0,0,0,0,0,0,0,0,0,0};
  string line;
  // skip the intro up until <init>
  do {
    getline(ifs,line);
  } while (line.compare("</init>") != 0);
  
  getline(ifs,line);
  
  int nevents=0,npass[4]={0,0,0,0};

  // get the event info
  while(getline(ifs,line)) { // A
    if(line.compare("<event>")==0) { // B
      nevents++;

      int idup, istup, mothup1, mothup2, icolup1, icolup2; 
      double  pupx, pupy, pupz, pupe, pupm, vtimup, spinup;
      
      TLorentzVector vq1,vq2,vq3,vq4;
      vq1.SetPxPyPzE  (0,0,0,0);
      vq2.SetPxPyPzE  (0,0,0,0);
      vq3.SetPxPyPzE  (0,0,0,0);
      vq4.SetPxPyPzE  (0,0,0,0);

      double pid1=0,pid2=0,pid3=0,pid4=0;

      // loop over the rest of them, skip initial line
      getline(ifs,line);
      getline(ifs,line);
      int pass[4] = {0,0,0,0}; int lType[3] = {0,0,0}; int nJets = 0;
      double weight = weightIni;
      while(line.compare("</event>") != 0) { // C
	stringstream sstmp(line);
	TString Line = line;
	if(!Line.Contains("#") && line.compare("") != 0) { // avoid crappy lines D
	  sstmp >> idup >> istup >> mothup1 >> mothup2 >> icolup1 >> icolup2 >> pupx >> pupy >> pupz >> pupe >> pupm >> vtimup >> spinup;
	  TLorentzVector vec;
	  vec.SetPxPyPzE(pupx,pupy,pupz,pupe);
	  
	  if(istup == 1){ // E
            if(TMath::Abs(idup) <6 || TMath::Abs(idup)==21) { // F
	      if     (vec.Pt() > vq1.Pt()){ //G
		vq4.SetPxPyPzE(vq3.Px(),vq3.Py(),vq3.Pz(),vq3.E());
		pid4=pid3;
		vq3.SetPxPyPzE(vq2.Px(),vq2.Py(),vq2.Pz(),vq2.E());
		pid3=pid2;
		vq2.SetPxPyPzE(vq1.Px(),vq1.Py(),vq1.Pz(),vq1.E());
		pid2=pid1;
		vq1.SetPxPyPzE(vec.Px(),vec.Py(),vec.Pz(),vec.E());
		pid1=idup;
	      } // G
	      else if(vec.Pt() > vq2.Pt()){ // G
		vq4.SetPxPyPzE(vq3.Px(),vq3.Py(),vq3.Pz(),vq3.E());
		pid4=pid3;
		vq3.SetPxPyPzE(vq2.Px(),vq2.Py(),vq2.Pz(),vq2.E());
		pid3=pid2;
		vq2.SetPxPyPzE(vec.Px(),vec.Py(),vec.Pz(),vec.E());
		pid2=idup;
	      } // G
	      else if(vec.Pt() > vq3.Pt()){ // G
		vq4.SetPxPyPzE(vq3.Px(),vq3.Py(),vq3.Pz(),vq3.E());
		pid4=pid3;
		vq3.SetPxPyPzE(vec.Px(),vec.Py(),vec.Pz(),vec.E());
		pid3=idup;
	      } // G
	      else if(vec.Pt() > vq4.Pt()){ // G 
		vq4.SetPxPyPzE(vec.Px(),vec.Py(),vec.Pz(),vec.E());
		pid4=idup;
	      } // G
            } // F
	    
	  } // E
	} // D
	getline(ifs,line);
      } // C
      //cout << endl;
      double ptq1,ptq2,ptq3,ptq4,etaq1,etaq2,etaq3,etaq4,phiq1,phiq2,phiq3,phiq4;
      if(vq1.Pt() > 0) { // C
	ptq1 = vq1.Pt();
	etaq1 = vq1.Eta();
	phiq1 = vq1.Phi();
      } 
      else            {
	ptq1 = 0.0;
	etaq1 = 0.0;
	phiq1 = 0.0;
      }
      if(vq2.Pt() > 0) {
	ptq2 = vq2.Pt();
	etaq2 = vq2.Eta();
	phiq2 = vq2.Phi();
      }
      else            {
	ptq2 = 0.0;
	etaq2 = 0.0;
	phiq2 = 0.0;
      }
      if(vq3.Pt() > 0) {
	ptq3 = vq3.Pt();
	etaq3 = vq3.Eta();
	phiq3 = vq3.Phi();
      }
      else            {
	ptq3 = 0.0;
	etaq3 = 0.0;
	phiq3 = 0.0;
      }
      if(vq4.Pt() > 0) {
	ptq4 = vq4.Pt();
	etaq4 = vq4.Eta();
	phiq4 = vq4.Phi();
      }
      else            {
	ptq4 = 0.0;
	etaq4 = 0.0;
	phiq4 = 0.0;
      } // C

      if (pid1==21) pid1=9;
      if (pid2==21) pid2=9;
      if (pid3==21) pid3=9;
      if (pid4==21) pid4=9;
      double pid=1000*fabs(pid1)+100*fabs(pid2)+10*fabs(pid3)+fabs(pid4);
      //cout << pid1 << " " << pid2 << " " << pid3 << " " << pid4 << " " << pid << endl;
      nt->Fill(ptq1,ptq2,ptq3,ptq4,etaq1,etaq2,etaq3,etaq4,phiq1,phiq2,phiq3,phiq4,pid,weight);
    } // B
    else { // B
      if(!TString(line).Contains("Random number generator exit values")) {
	//ofs << line << endl;
      }
    } // B
  } // A
  
  ifs.close();
  nt->Write();
  outtuple->Close();

}

double DeltaPhi(double phi1, double phi2)
{
  // Compute DeltaPhi between two given angles. Results is in [-pi/2,pi/2].
  double dphi = TMath::Abs(phi1-phi2);
  while (dphi>TMath::Pi())
    dphi = TMath::Abs(dphi - TMath::TwoPi());
  return(dphi);
}
