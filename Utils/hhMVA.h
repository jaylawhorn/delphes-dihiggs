#ifndef HHMVA_H
#define HHMVA_H

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TString.h>
#include <TMath.h>
#include <string>
#include <vector>
#include <utility>
#include <map>
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

using namespace std;
using namespace TMVA;

class hhMVA {
 public:
  hhMVA();
  ~hhMVA();
  
  enum MVAType {
    kTauTau=1,
    kMuTau=2,
    kElTau=3,
    kElMu=4,
  };
  
  void Intialize(MVAType type);
  
  Float_t GetBDTValue(Float_t ptTau1, Float_t ptTau2, Float_t ptB1, Float_t ptB2, 
		      Float_t mTT, Float_t ptTT, Float_t mBB, Float_t ptBB,
		      Float_t mHH, Float_t ptHH, Float_t mt2pileup,
		      Float_t dRbb, Float_t dRtt, Float_t dRhh);
  
 protected:
  TMVA::Reader         *fTMVAReader;
  MVAType              fMVAType;
  
  Float_t             fMVAVar_ptTau1;
  Float_t             fMVAVar_ptTau2;
  Float_t             fMVAVar_ptB1;
  Float_t             fMVAVar_ptB2;
  Float_t             fMVAVar_mTT;
  Float_t             fMVAVar_ptTT;
  Float_t             fMVAVar_mBB1;
  Float_t             fMVAVar_ptBB1;
  Float_t             fMVAVar_mHH;
  Float_t             fMVAVar_ptHH;
  Float_t             fMVAVar_mt2pileup;
  Float_t             fMVAVar_dRbb;
  Float_t             fMVAVar_dRtt;
  Float_t             fMVAVar_dRhh;

  TString              fTTWeights;
  TString              fMTWeights;
  TString              fETWeights;
  TString              fEMWeights;
  
};
#endif
  
