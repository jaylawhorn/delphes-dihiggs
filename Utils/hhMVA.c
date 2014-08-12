#include "hhMVA.h"

#include "TMVA/Reader.h"
#include "TMVA/Tools.h"

using namespace std;
using namespace TMVA;

hhMVA::hhMVA() {
  fTMVAReader = 0;
  //fTTWeights = "../MVA/TMVA_BDT_tt.weights.xml";
  fMTWeights = "../MVA/TMVAClassification_BDT_mt.weights.xml";
  fETWeights = "../MVA/TMVAClassification_BDT_et.weights.xml";
  fEMWeights = "../MVA/TMVAClassification_BDT_em.weights.xml";
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

  if (fMVAType==kMuTau || fMVAType==kElTau) { 
    fTMVAReader->AddVariable( "ptTau1",     &fMVAVar_ptTau1   );
    fTMVAReader->AddVariable( "ptTau2",     &fMVAVar_ptTau2   );
    fTMVAReader->AddVariable( "ptB1",       &fMVAVar_ptB1     );
    fTMVAReader->AddVariable( "ptB2",       &fMVAVar_ptB2     );
    fTMVAReader->AddVariable( "mTT",        &fMVAVar_mTT      );
    fTMVAReader->AddVariable( "ptTT",       &fMVAVar_ptTT     );
    fTMVAReader->AddVariable( "mBB1",       &fMVAVar_mBB1     );
    fTMVAReader->AddVariable( "ptBB1",      &fMVAVar_ptBB1    );
    fTMVAReader->AddVariable( "mHH",        &fMVAVar_mHH      );
    fTMVAReader->AddVariable( "ptHH",       &fMVAVar_ptHH     );
    fTMVAReader->AddVariable( "mt2pileup",  &fMVAVar_mt2pileup);
    fTMVAReader->AddVariable( "dRbb := sqrt((etaB1-etaB2)**2+(phiB1-phiB2)**2)",         &fMVAVar_dRbb );
    fTMVAReader->AddVariable( "dRtt := sqrt((etaTau1-etaTau2)**2+(phiTau1-phiTau2)**2)", &fMVAVar_dRtt );
    fTMVAReader->AddVariable( "dRhh := sqrt((etaBB1-etaTT)**2+(phiBB1-phiTT)**2)",       &fMVAVar_dRhh );
  }
  else if (fMVAType==kElMu) {
    fTMVAReader->AddVariable( "mTT",        &fMVAVar_mTT      );
    fTMVAReader->AddVariable( "mBB1",       &fMVAVar_mBB1     );
    fTMVAReader->AddVariable( "mt2pileup",  &fMVAVar_mt2pileup);
    fTMVAReader->AddVariable( "dRbb := sqrt((etaB1-etaB2)**2+(phiB1-phiB2)**2)",         &fMVAVar_dRbb );
    fTMVAReader->AddVariable( "dRtt := sqrt((etaTau1-etaTau2)**2+(phiTau1-phiTau2)**2)", &fMVAVar_dRtt );
    fTMVAReader->AddVariable( "dRhh := sqrt((etaBB1-etaTT)**2+(phiBB1-phiTT)**2)",       &fMVAVar_dRhh );
  }
  
  //if (fMVAType==kTauTau) weightfile=fTTWeights;
  if (fMVAType==kMuTau)  weightfile=fMTWeights;
  if (fMVAType==kElTau)  weightfile=fETWeights;
  if (fMVAType==kElMu)   weightfile=fEMWeights;

  fTMVAReader->BookMVA("BDT method", weightfile);

  assert(fTMVAReader);

}

Float_t hhMVA::GetBDTValue(Float_t ptTau1, Float_t ptTau2, Float_t ptB1, Float_t ptB2,
			   Float_t mTT, Float_t ptTT, Float_t mBB, Float_t ptBB,
			   Float_t mHH, Float_t ptHH, Float_t mt2pileup,
			   Float_t dRbb, Float_t dRtt, Float_t dRhh) {
  
  if (!fTMVAReader) {
    cout << "TMVA reader not initialized properly" << endl;
    return -999;
  }

  if (fMVAType==kMuTau || fMVAType==kElTau) {
    fMVAVar_ptTau1=ptTau1;
    fMVAVar_ptTau2=ptTau2;
    fMVAVar_ptB1  =ptB1;
    fMVAVar_ptB2  =ptB2;
    fMVAVar_mTT   =mTT;
    fMVAVar_ptTT  =ptTT;
    fMVAVar_mBB1  =mBB;
    fMVAVar_ptBB1 =ptBB;
    fMVAVar_mHH   =mHH;
    fMVAVar_ptHH  =ptHH;
    fMVAVar_mt2pileup=mt2pileup;
    fMVAVar_dRbb  =dRbb;
    fMVAVar_dRtt  =dRtt;
    fMVAVar_dRhh  =dRhh;
  }
  else if (fMVAType==kElMu) {
    fMVAVar_mTT   =mTT;
    fMVAVar_mBB1  =mBB;
    fMVAVar_mt2pileup=mt2pileup;
    fMVAVar_dRbb  =dRbb;
    fMVAVar_dRtt  =dRtt;
    fMVAVar_dRhh  =dRhh;
  }
  
  TMVA::Reader *reader = 0;
  reader = fTMVAReader;
  
  return reader->EvaluateMVA("BDT method");
  
}
