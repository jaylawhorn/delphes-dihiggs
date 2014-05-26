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

void getall(const TString inputfile) { 

  Int_t totalEvents=0;
  Int_t nEvents=0;

  ifstream ifs;
  ifs.open(inputfile.Data()); assert(ifs.is_open());
  string line;
  while(getline(ifs,line)) {
    string fname;
    stringstream ss(line);
    ss >> fname;

    TChain chain("Delphes");
    chain.Add(TString(fname));
    ExRootTreeReader treeReader(&chain);
    nEvents=treeReader.GetEntries();
    cout << nEvents << endl;
    totalEvents+=nEvents;

  }
  ifs.close();

  cout << totalEvents << endl;

}
