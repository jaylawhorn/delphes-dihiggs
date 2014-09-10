#include "hhMVA.h"

#include "TMVA/Reader.h"
#include "TMVA/Tools.h"

using namespace std;
using namespace TMVA;

hhMVA::hhMVA() {
  fTMVAReader = 0;
  fTTWeights = "../MVA/TMVAClassification_BDT.tt.weights.xml";
  fMTWeights = "../MVA/TMVAClassification_BDT.mt.weights.xml";
  fETWeights = "../MVA/TMVAClassification_BDT.et.weights.xml";
  fEMWeights = "../MVA/TMVAClassification_BDT.em.weights.xml";
}

hhMVA::~hhMVA() {
  if (fTMVAReader) delete fTMVAReader;
}

void hhMVA::Intialize(MVAType type) {
  fMVAType=type;

  if (fTMVAReader) delete fTMVAReader;

  fTMVAReader = new TMVA::Reader( "!Color:!Silent:Error" );
  fTMVAReader->SetVerbose(kTRUE);

  TString weightfile;

  fTMVAReader->AddVariable( "mTT",        &fMVAVar_mTT      );
  fTMVAReader->AddVariable( "ptTT",       &fMVAVar_ptTT     );
  fTMVAReader->AddVariable( "mBB1",       &fMVAVar_mBB1     );
  fTMVAReader->AddVariable( "ptBB1",      &fMVAVar_ptBB1    );
  fTMVAReader->AddVariable( "mHH",        &fMVAVar_mHH      );
  fTMVAReader->AddVariable( "ptHH",       &fMVAVar_ptHH     );
  fTMVAReader->AddVariable( "mt2pileup",  &fMVAVar_mt2      );
  
  fTMVAReader->AddVariable( "dRBB1",      &fMVAVar_dRbb     );
  fTMVAReader->AddVariable( "dRTT",       &fMVAVar_dRtt     );
  fTMVAReader->AddVariable( "dRHH",       &fMVAVar_dRhh     );
  
  if (fMVAType==kTauTau) weightfile=fTTWeights;
  if (fMVAType==kMuTau)  weightfile=fMTWeights;
  if (fMVAType==kElTau)  weightfile=fETWeights;
  if (fMVAType==kElMu)   weightfile=fEMWeights;
  
  fTMVAReader->BookMVA("BDT method", weightfile);

  assert(fTMVAReader);

}

Float_t hhMVA::GetBDTValue(Float_t mTT, Float_t ptTT, Float_t mBB, Float_t ptBB,
			   Float_t mHH, Float_t ptHH, Float_t mt2,
			   Float_t dRbb, Float_t dRtt, Float_t dRhh) {
  
  if (!fTMVAReader) {
    cout << "TMVA reader not initialized properly" << endl;
    return -999;
  }

  fMVAVar_mTT   =mTT;
  fMVAVar_ptTT  =ptTT;
  fMVAVar_mBB1  =mBB;
  fMVAVar_ptBB1 =ptBB;
  fMVAVar_mHH   =mHH;
  fMVAVar_ptHH  =ptHH;
  fMVAVar_mt2   =mt2;
  fMVAVar_dRbb  =dRbb;
  fMVAVar_dRtt  =dRtt;
  fMVAVar_dRhh  =dRhh;

  TMVA::Reader *reader = 0;
  reader = fTMVAReader;
  
  return reader->EvaluateMVA("BDT method");
  
}
