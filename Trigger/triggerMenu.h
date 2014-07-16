#ifndef TRIGGER_MENU_H
#define TRIGGER_MENU_H

namespace Trigger {
  enum {
    kSingleMu=1, 
    kDoubleMu=2,
    kEleMu=4,
    kSingleEle=8,
    kSingleIsoEle=16,
    kSingleGamma=32, 
    kIsoEleEle=64,
    kDoubleGamma=128,
    kSingleTau=256,
    kDoubleTau=512, 
    kEleTau=1024,
    kTauMu=2048,
    kSingleJet=4096,
    kDoubleJet=8192,
    kQuadJet=16384,
    kEleJet=32768,
    kMuJet=65536,
    kEleMet=131072,
    kMuMet=262144,
    kHt=524288
  };
  
  struct Event {
    UInt_t isReco;
    Float_t e1pt, e1eta, e1phi;
    Float_t e2pt, e2eta, e2phi;
    Float_t m1pt, m1eta, m1phi;
    Float_t m2pt, m2eta, m2phi;
    Float_t g1pt, g1eta, g1phi;
    Float_t g2pt, g2eta, g2phi;
    Float_t t1pt, t1eta, t1phi;
    Float_t t2pt, t2eta, t2phi;
    Float_t j1pt, j1eta, j1phi;
    Float_t j2pt, j2eta, j2phi;
    Float_t j3pt, j3eta, j3phi;
    Float_t j4pt, j4eta, j4phi;
    Float_t mht, ht;
    UInt_t triggerBits;
  };
  
  UInt_t checkTriggers(Trigger::Event data) {
    UInt_t result=0;
    if ( (data.m1pt>21.1 && abs(data.m1eta)<2.4) || (data.m2pt>21.1 && abs(data.m2eta)<2.4) )
      result+=Trigger::kSingleMu;
    
    if ( data.m1pt>18.4 && abs(data.m1eta)<2.4 && data.m2pt>13.0 && abs(data.m2eta)<2.4)
      result+=Trigger::kDoubleMu;
    
    if ( ( (data.e1pt>19.6 && abs(data.e1eta)<2.4) || (data.e2pt>19.6 && abs(data.e2eta)<2.4) ) && ( (data.m1pt>10.9 && abs(data.m1eta)<2.4) || (data.m2pt>10.9 && abs(data.m2eta)<2.4) ) )
      result+=Trigger::kEleMu;
    
    if ( (data.e1pt>35.3 && abs(data.e1eta)<2.4) || (data.e2pt>35.3 && abs(data.e1eta)<2.4) )
      result+=Trigger::kSingleEle;

    if ( (data.e1pt>31.2 && abs(data.e1eta)<2.4) || (data.e2pt>31.2 && abs(data.e1eta)<2.4) )
      // need to implement iso check at gen level
      result+=Trigger::kSingleIsoEle;
    
    if ( (data.g1pt>34.8 && abs(data.g1eta)<2.4) || (data.g2pt>34.8 && abs(data.g2eta)<2.4) )
      result+=Trigger::kSingleGamma;
    
    if ( (data.e1pt>23.3 && abs(data.e1eta)<2.4) && ( (data.e2pt>16.7 && abs(data.e2eta)<2.4) || (data.g1pt>16.7 && abs(data.g1eta)<2.4) || (data.g2pt>16.7 && abs(data.g2eta)<2.4) ) )
      result+=Trigger::kIsoEleEle;
    
    if ( data.g1pt>23.3 && abs(data.g1eta)<2.4 && data.g2pt>14.7 && abs(data.g2eta)<2.4)
      result+=Trigger::kDoubleGamma;

    if ( (data.t1pt>101 && abs(data.t1eta)<2.4) || (data.t2pt>101 && abs(data.t2eta)<2.4) )
      result+=Trigger::kSingleTau;
    
    if ( data.t1pt>63.8 && abs(data.t1eta)<2.4 && data.t2pt>62.6 && abs(data.t2eta)<2.4 )
      result+=Trigger::kDoubleTau;
    
    if ( (data.e1pt>20.4 && abs(data.e1eta)<2.4) || (data.e2pt>20.4 && abs(data.e2eta)<2.4) ) {
      if ( (data.t1pt>53.6 && abs(data.t1eta)<2.4) || (data.t1pt>53.6 && abs(data.t1eta)<2.4) ) {
	result+=Trigger::kEleTau;
      }
    }
    
    if ( (data.m1pt>16.0 && abs(data.m1eta)<2.4) || (data.m2pt>16.0 && abs(data.m2eta)<2.4) ) {
      if ( (data.t1pt>51.9 && abs(data.t1eta)<2.4) || (data.t1pt>51.9 && abs(data.t1eta)<2.4) ) {
	result+=Trigger::kTauMu;
      }
    }
    
    return result;
  }
  
}

Float_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 ) {
  const Float_t pi = 3.14159265358979;

  Float_t etaDiff = (eta1-eta2);
  Float_t phiDiff = fabs(phi1-phi2);
  while ( phiDiff>pi ) phiDiff = fabs(phiDiff-2.0*pi);

  Float_t deltaRSquared = etaDiff*etaDiff + phiDiff*phiDiff;

  return TMath::Sqrt(deltaRSquared);

}


#endif
