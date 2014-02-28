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
#include "bJetScaleCorr.h"
#endif

void testJetCorr() {


  TF1 test2 = (*f(0.15));
  cout << test2.GetName() << endl;
  TF1 test3 = (*f(0.45));
  cout << test3.GetName() << endl;
  TF1 test4 = (*f(0.75));
  cout << test4.GetName() << endl;
  TF1 test5 = (*f(1.05));
  cout << test5.GetName() << endl;
  TF1 test6 = (*f(1.35));
  cout << test6.GetName() << endl;
  TF1 test7 = (*f(1.65));
  cout << test7.GetName() << endl;
  TF1 test8 = (*f(1.95));
  cout << test8.GetName() << endl;
  TF1 test9 = (*f(2.25));
  cout << test9.GetName() << endl;
  TF1 test10 = (*f(2.55));
  cout << test10.GetName() << endl;
  TF1 test11 = (*f(2.85));
  cout << test11.GetName() << endl;
  TF1 test12 = (*f(3.5));
  cout << test12.GetName() << endl;
  TF1 test13 = (*f(4.5));
  cout << test13.GetName() << endl;

}
