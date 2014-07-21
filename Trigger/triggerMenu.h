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
    UInt_t triggerBits126;
    UInt_t triggerBits180;
    UInt_t triggerBits250;
    UInt_t triggerBits350;
  };

  // from RateTables_v4
  UInt_t checkTriggers126(Trigger::Event data) {
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

    if ( (data.j1pt>188 && abs(data.j1eta)<2.4) || (data.j2pt>188 && abs(data.j2eta)<2.4) || (data.j3pt>188 && abs(data.j3eta)<2.4) || (data.j4pt>188 && abs(data.j4eta)<2.4) ) 
      result+=Trigger::kSingleJet;

    if (data.j1pt>138 && abs(data.j1eta)<2.4) {
      if ( (data.j2pt>138 && abs(data.j2eta)<2.4) || (data.j3pt>138 && abs(data.j3eta)<2.4) || (data.j4pt>138 && abs(data.j4eta)<2.4) ) 
	result+=Trigger::kDoubleJet;
    }
    else if (data.j2pt>138 && abs(data.j2eta)<2.4) {
      if ( (data.j3pt>138 && abs(data.j3eta)<2.4) || (data.j4pt>138 && abs(data.j4eta)<2.4) )
        result+=Trigger::kDoubleJet;
    }
    else if (data.j3pt>138 && abs(data.j3eta)<2.4) {
      if (data.j4pt>138 && abs(data.j4eta)<2.4)
        result+=Trigger::kDoubleJet;
    }

    if ( (data.j1pt>76.5 && abs(data.j1eta)<2.4) && (data.j2pt>76.5 && abs(data.j2eta)<2.4) && (data.j3pt>76.5 && abs(data.j3eta)<2.4) && (data.j4pt>76.5 && abs(data.j4eta)<2.4) )
      result+=Trigger::kQuadJet;

    if ( (data.e1pt>26.3 && abs(data.e1eta)<2.4) || (data.e2pt>26.3 && abs(data.e2eta)<2.4) ) {
      if ( (data.j1pt>73.1 && abs(data.j1eta)<2.4) || (data.j2pt>73.1 && abs(data.j2eta)<2.4) || (data.j3pt>73.1 && abs(data.j3eta)<2.4) || (data.j4pt>73.1 && abs(data.j4eta)<2.4) )
	result+=Trigger::kEleJet;
    }

    if ( (data.m1pt>18.1 && abs(data.m1eta)<2.4) || (data.m2pt>18.1 && abs(data.m2eta)<2.4) ) {
      if ( (data.j1pt>71.9 && abs(data.j1eta)<2.4) || (data.j2pt>71.9 && abs(data.j2eta)<2.4) || (data.j3pt>71.9 && abs(data.j3eta)<2.4) || (data.j4pt>71.9 && abs(data.j4eta)<2.4) )
	result+=Trigger::kMuJet;
    }

    if ( (data.e1pt>24.8 && abs(data.e1eta)<2.4) || (data.e2pt>24.8 && abs(data.e2eta)<2.4) ) {
      if (data.mht>102) {
	result+=Trigger::kEleMet;
      }
    }
    
    if ( (data.m1pt>15.7 && abs(data.m1eta)<2.4) || (data.m2pt>15.7 && abs(data.m2eta)<2.4) ) {
      if (data.mht>68.6) {
	result+=Trigger::kMuMet;
      }
    }
    
    if (data.ht>384)
      result+=Trigger::kHt;

    return result;
  }

  // from approval talk
  UInt_t checkTriggers180(Trigger::Event data) {
    UInt_t result=0;
    if ( (data.m1pt>18 && abs(data.m1eta)<2.4) || (data.m2pt>18 && abs(data.m2eta)<2.4) )
      result+=Trigger::kSingleMu;
    
    if ( data.m1pt>14 && abs(data.m1eta)<2.4 && data.m2pt>10 && abs(data.m2eta)<2.4)
      result+=Trigger::kDoubleMu;
    
    if ( ( (data.e1pt>19 && abs(data.e1eta)<2.4) || (data.e2pt>19 && abs(data.e2eta)<2.4) ) && ( (data.m1pt>10.5 && abs(data.m1eta)<2.4) || (data.m2pt>10.5 && abs(data.m2eta)<2.4) ) )
      result+=Trigger::kEleMu;
    
    if ( (data.e1pt>31 && abs(data.e1eta)<2.4) || (data.e2pt>31 && abs(data.e1eta)<2.4) )
      result+=Trigger::kSingleEle;

    if ( (data.e1pt>27 && abs(data.e1eta)<2.4) || (data.e2pt>27 && abs(data.e1eta)<2.4) )
      result+=Trigger::kSingleIsoEle;
    
    if ( (data.g1pt>31 && abs(data.g1eta)<2.4) || (data.g2pt>31 && abs(data.g2eta)<2.4) )
      result+=Trigger::kSingleGamma;
    
    if ( (data.e1pt>22 && abs(data.e1eta)<2.4) && ( (data.e2pt>16 && abs(data.e2eta)<2.4) || (data.g1pt>16 && abs(data.g1eta)<2.4) || (data.g2pt>16 && abs(data.g2eta)<2.4) ) )
      result+=Trigger::kIsoEleEle;
    
    if ( data.g1pt>22 && abs(data.g1eta)<2.4 && data.g2pt>16 && abs(data.g2eta)<2.4)
      result+=Trigger::kDoubleGamma;

    if ( (data.t1pt>88 && abs(data.t1eta)<2.4) || (data.t2pt>88 && abs(data.t2eta)<2.4) )
      result+=Trigger::kSingleTau;
    
    if ( data.t1pt>56 && abs(data.t1eta)<2.4 && data.t2pt>56 && abs(data.t2eta)<2.4 )
      result+=Trigger::kDoubleTau;
    
    if ( (data.e1pt>19 && abs(data.e1eta)<2.4) || (data.e2pt>19 && abs(data.e2eta)<2.4) ) {
      if ( (data.t1pt>50 && abs(data.t1eta)<2.4) || (data.t1pt>50 && abs(data.t1eta)<2.4) ) {
	result+=Trigger::kEleTau;
      }
    }
    
    if ( (data.m1pt>14 && abs(data.m1eta)<2.4) || (data.m2pt>14 && abs(data.m2eta)<2.4) ) {
      if ( (data.t1pt>45 && abs(data.t1eta)<2.4) || (data.t1pt>45 && abs(data.t1eta)<2.4) ) {
	result+=Trigger::kTauMu;
      }
    }

    if ( (data.j1pt>173 && abs(data.j1eta)<2.4) || (data.j2pt>173 && abs(data.j2eta)<2.4) || (data.j3pt>173 && abs(data.j3eta)<2.4) || (data.j4pt>173 && abs(data.j4eta)<2.4) ) 
      result+=Trigger::kSingleJet;

    if (data.j1pt>125 && abs(data.j1eta)<2.4) {
      if ( (data.j2pt>125 && abs(data.j2eta)<2.4) || (data.j3pt>125 && abs(data.j3eta)<2.4) || (data.j4pt>125 && abs(data.j4eta)<2.4) ) 
	result+=Trigger::kDoubleJet;
    }
    else if (data.j2pt>125 && abs(data.j2eta)<2.4) {
      if ( (data.j3pt>125 && abs(data.j3eta)<2.4) || (data.j4pt>125 && abs(data.j4eta)<2.4) )
        result+=Trigger::kDoubleJet;
    }
    else if (data.j3pt>125 && abs(data.j3eta)<2.4) {
      if (data.j4pt>125 && abs(data.j4eta)<2.4)
        result+=Trigger::kDoubleJet;
    }

    if ( (data.j1pt>72 && abs(data.j1eta)<2.4) && (data.j2pt>72 && abs(data.j2eta)<2.4) && (data.j3pt>72 && abs(data.j3eta)<2.4) && (data.j4pt>72 && abs(data.j4eta)<2.4) )
      result+=Trigger::kQuadJet;

    if ( (data.e1pt>23 && abs(data.e1eta)<2.4) || (data.e2pt>23 && abs(data.e2eta)<2.4) ) {
      if ( (data.j1pt>66 && abs(data.j1eta)<2.4) || (data.j2pt>66 && abs(data.j2eta)<2.4) || (data.j3pt>66 && abs(data.j3eta)<2.4) || (data.j4pt>66 && abs(data.j4eta)<2.4) )
	result+=Trigger::kEleJet;
    }

    if ( (data.m1pt>16 && abs(data.m1eta)<2.4) || (data.m2pt>16 && abs(data.m2eta)<2.4) ) {
      if ( (data.j1pt>66 && abs(data.j1eta)<2.4) || (data.j2pt>66 && abs(data.j2eta)<2.4) || (data.j3pt>66 && abs(data.j3eta)<2.4) || (data.j4pt>66 && abs(data.j4eta)<2.4) )
	result+=Trigger::kMuJet;
    }

    if ( (data.e1pt>23 && abs(data.e1eta)<2.4) || (data.e2pt>23 && abs(data.e2eta)<2.4) ) {
      if (data.mht>95) {
	result+=Trigger::kEleMet;
      }
    }

    if ( (data.m1pt>16 && abs(data.m1eta)<2.4) || (data.m2pt>16 && abs(data.m2eta)<2.4) ) {
      if (data.mht>95) {
	result+=Trigger::kMuMet;
      }
    }

    if (data.ht>350)
      result+=Trigger::kHt;
    
    return result;
  }

  // from RateTables_v2
  UInt_t checkTriggers250(Trigger::Event data) {
    UInt_t result=0;
    if ( (data.m1pt>17.3 && abs(data.m1eta)<2.4) || (data.m2pt>17.3 && abs(data.m2eta)<2.4) )
      result+=Trigger::kSingleMu;
    
    if ( data.m1pt>11.9 && abs(data.m1eta)<2.4 && data.m2pt>8.6 && abs(data.m2eta)<2.4)
      result+=Trigger::kDoubleMu;
    
    if ( ( (data.e1pt>18.5 && abs(data.e1eta)<2.4) || (data.e2pt>18.5 && abs(data.e2eta)<2.4) ) && ( (data.m1pt>10 && abs(data.m1eta)<2.4) || (data.m2pt>10 && abs(data.m2eta)<2.4) ) )
      result+=Trigger::kEleMu;
    
    if ( (data.e1pt>29.6 && abs(data.e1eta)<2.4) || (data.e2pt>29.6 && abs(data.e1eta)<2.4) )
      result+=Trigger::kSingleEle;

    if ( (data.e1pt>25.1 && abs(data.e1eta)<2.4) || (data.e2pt>25.1 && abs(data.e1eta)<2.4) )
      result+=Trigger::kSingleIsoEle;
    
    if ( (data.g1pt>29.6 && abs(data.g1eta)<2.4) || (data.g2pt>29.6 && abs(data.g2eta)<2.4) )
      result+=Trigger::kSingleGamma;
    
    if ( (data.e1pt>21.5 && abs(data.e1eta)<2.4) && ( (data.e2pt>15 && abs(data.e2eta)<2.4) || (data.g1pt>15 && abs(data.g1eta)<2.4) || (data.g2pt>15 && abs(data.g2eta)<2.4) ) )
      result+=Trigger::kIsoEleEle;
    
    if ( data.g1pt>21.7 && abs(data.g1eta)<2.4 && data.g2pt>13 && abs(data.g2eta)<2.4)
      result+=Trigger::kDoubleGamma;

    if ( (data.t1pt>83.9 && abs(data.t1eta)<2.4) || (data.t2pt>83.9 && abs(data.t2eta)<2.4) )
      result+=Trigger::kSingleTau;
    
    if ( data.t1pt>54 && abs(data.t1eta)<2.4 && data.t2pt>54 && abs(data.t2eta)<2.4 )
      result+=Trigger::kDoubleTau;
    
    if ( (data.e1pt>18.5 && abs(data.e1eta)<2.4) || (data.e2pt>18.5 && abs(data.e2eta)<2.4) ) {
      if ( (data.t1pt>48 && abs(data.t1eta)<2.4) || (data.t1pt>48 && abs(data.t1eta)<2.4) ) {
	result+=Trigger::kEleTau;
      }
    }
    
    if ( (data.m1pt>13 && abs(data.m1eta)<2.4) || (data.m2pt>13 && abs(data.m2eta)<2.4) ) {
      if ( (data.t1pt>43.4 && abs(data.t1eta)<2.4) || (data.t1pt>43.4 && abs(data.t1eta)<2.4) ) {
	result+=Trigger::kTauMu;
      }
    }

    if ( (data.j1pt>167 && abs(data.j1eta)<2.4) || (data.j2pt>167 && abs(data.j2eta)<2.4) || (data.j3pt>167 && abs(data.j3eta)<2.4) || (data.j4pt>167 && abs(data.j4eta)<2.4) ) 
      result+=Trigger::kSingleJet;

    if (data.j1pt>121 && abs(data.j1eta)<2.4) {
      if ( (data.j2pt>121 && abs(data.j2eta)<2.4) || (data.j3pt>121 && abs(data.j3eta)<2.4) || (data.j4pt>121 && abs(data.j4eta)<2.4) ) 
	result+=Trigger::kDoubleJet;
    }
    else if (data.j2pt>121 && abs(data.j2eta)<2.4) {
      if ( (data.j3pt>121 && abs(data.j3eta)<2.4) || (data.j4pt>121 && abs(data.j4eta)<2.4) )
        result+=Trigger::kDoubleJet;
    }
    else if (data.j3pt>121 && abs(data.j3eta)<2.4) {
      if (data.j4pt>121 && abs(data.j4eta)<2.4)
        result+=Trigger::kDoubleJet;
    }

    if ( (data.j1pt>71 && abs(data.j1eta)<2.4) && (data.j2pt>71 && abs(data.j2eta)<2.4) && (data.j3pt>71 && abs(data.j3eta)<2.4) && (data.j4pt>71 && abs(data.j4eta)<2.4) )
      result+=Trigger::kQuadJet;

    if ( (data.e1pt>22 && abs(data.e1eta)<2.4) || (data.e2pt>22 && abs(data.e2eta)<2.4) ) {
      if ( (data.j1pt>64 && abs(data.j1eta)<2.4) || (data.j2pt>64 && abs(data.j2eta)<2.4) || (data.j3pt>64 && abs(data.j3eta)<2.4) || (data.j4pt>64 && abs(data.j4eta)<2.4) )
	result+=Trigger::kEleJet;
    }

    if ( (data.m1pt>15.4 && abs(data.m1eta)<2.4) || (data.m2pt>15.4 && abs(data.m2eta)<2.4) ) {
      if ( (data.j1pt>64 && abs(data.j1eta)<2.4) || (data.j2pt>64 && abs(data.j2eta)<2.4) || (data.j3pt>64 && abs(data.j3eta)<2.4) || (data.j4pt>64 && abs(data.j4eta)<2.4) )
	result+=Trigger::kMuJet;
    }

    if ( (data.e1pt>22.3 && abs(data.e1eta)<2.4) || (data.e2pt>22.3 && abs(data.e2eta)<2.4) ) {
      if (data.mht>92) {
	result+=Trigger::kEleMet;
      }
    }

    if ( (data.m1pt>13.9 && abs(data.m1eta)<2.4) || (data.m2pt>13.9 && abs(data.m2eta)<2.4) ) {
      if (data.mht>61) {
	result+=Trigger::kMuMet;
      }
    }

    if (data.ht>340)
      result+=Trigger::kHt;
    
    return result;
  }

  // from RateTables_v2
  UInt_t checkTriggers350(Trigger::Event data) {
    UInt_t result=0;
    if ( (data.m1pt>15.5 && abs(data.m1eta)<2.4) || (data.m2pt>15.5 && abs(data.m2eta)<2.4) )
      result+=Trigger::kSingleMu;
    
    if ( data.m1pt>10.4 && abs(data.m1eta)<2.4 && data.m2pt>7.6 && abs(data.m2eta)<2.4)
      result+=Trigger::kDoubleMu;
    
    if ( ( (data.e1pt>16.3 && abs(data.e1eta)<2.4) || (data.e2pt>16.3 && abs(data.e2eta)<2.4) ) && ( (data.m1pt>8.8 && abs(data.m1eta)<2.4) || (data.m2pt>8.8 && abs(data.m2eta)<2.4) ) )
      result+=Trigger::kEleMu;
    
    if ( (data.e1pt>26.3 && abs(data.e1eta)<2.4) || (data.e2pt>26.3 && abs(data.e1eta)<2.4) )
      result+=Trigger::kSingleEle;

    if ( (data.e1pt>22.8 && abs(data.e1eta)<2.4) || (data.e2pt>22.8 && abs(data.e1eta)<2.4) )
      result+=Trigger::kSingleIsoEle;
    
    if ( (data.g1pt>27.3 && abs(data.g1eta)<2.4) || (data.g2pt>27.3 && abs(data.g2eta)<2.4) )
      result+=Trigger::kSingleGamma;
    
    if ( (data.e1pt>20.6 && abs(data.e1eta)<2.4) && ( (data.e2pt>15 && abs(data.e2eta)<2.4) || (data.g1pt>15 && abs(data.g1eta)<2.4) || (data.g2pt>15 && abs(data.g2eta)<2.4) ) )
      result+=Trigger::kIsoEleEle;
    
    if ( data.g1pt>20.7 && abs(data.g1eta)<2.4 && data.g2pt>13 && abs(data.g2eta)<2.4)
      result+=Trigger::kDoubleGamma;

    if ( (data.t1pt>76 && abs(data.t1eta)<2.4) || (data.t2pt>76 && abs(data.t2eta)<2.4) )
      result+=Trigger::kSingleTau;
    
    if ( data.t1pt>50 && abs(data.t1eta)<2.4 && data.t2pt>50 && abs(data.t2eta)<2.4 )
      result+=Trigger::kDoubleTau;
    
    if ( (data.e1pt>17.7 && abs(data.e1eta)<2.4) || (data.e2pt>17.7 && abs(data.e2eta)<2.4) ) {
      if ( (data.t1pt>46 && abs(data.t1eta)<2.4) || (data.t1pt>46 && abs(data.t1eta)<2.4) ) {
	result+=Trigger::kEleTau;
      }
    }
    
    if ( (data.m1pt>12 && abs(data.m1eta)<2.4) || (data.m2pt>12 && abs(data.m2eta)<2.4) ) {
      if ( (data.t1pt>40.8 && abs(data.t1eta)<2.4) || (data.t1pt>40.8 && abs(data.t1eta)<2.4) ) {
	result+=Trigger::kTauMu;
      }
    }

    if ( (data.j1pt>158 && abs(data.j1eta)<2.4) || (data.j2pt>158 && abs(data.j2eta)<2.4) || (data.j3pt>158 && abs(data.j3eta)<2.4) || (data.j4pt>158 && abs(data.j4eta)<2.4) ) 
      result+=Trigger::kSingleJet;

    if (data.j1pt>114 && abs(data.j1eta)<2.4) {
      if ( (data.j2pt>114 && abs(data.j2eta)<2.4) || (data.j3pt>114 && abs(data.j3eta)<2.4) || (data.j4pt>114 && abs(data.j4eta)<2.4) ) 
	result+=Trigger::kDoubleJet;
    }
    else if (data.j2pt>114 && abs(data.j2eta)<2.4) {
      if ( (data.j3pt>114 && abs(data.j3eta)<2.4) || (data.j4pt>114 && abs(data.j4eta)<2.4) )
        result+=Trigger::kDoubleJet;
    }
    else if (data.j3pt>114 && abs(data.j3eta)<2.4) {
      if (data.j4pt>114 && abs(data.j4eta)<2.4)
        result+=Trigger::kDoubleJet;
    }

    if ( (data.j1pt>69 && abs(data.j1eta)<2.4) && (data.j2pt>69 && abs(data.j2eta)<2.4) && (data.j3pt>69 && abs(data.j3eta)<2.4) && (data.j4pt>69 && abs(data.j4eta)<2.4) )
      result+=Trigger::kQuadJet;

    if ( (data.e1pt>20.4 && abs(data.e1eta)<2.4) || (data.e2pt>20.4 && abs(data.e2eta)<2.4) ) {
      if ( (data.j1pt>60 && abs(data.j1eta)<2.4) || (data.j2pt>60 && abs(data.j2eta)<2.4) || (data.j3pt>60 && abs(data.j3eta)<2.4) || (data.j4pt>60 && abs(data.j4eta)<2.4) )
	result+=Trigger::kEleJet;
    }

    if ( (data.m1pt>14.3 && abs(data.m1eta)<2.4) || (data.m2pt>14.3 && abs(data.m2eta)<2.4) ) {
      if ( (data.j1pt>61 && abs(data.j1eta)<2.4) || (data.j2pt>61 && abs(data.j2eta)<2.4) || (data.j3pt>61 && abs(data.j3eta)<2.4) || (data.j4pt>61 && abs(data.j4eta)<2.4) )
	result+=Trigger::kMuJet;
    }

    if ( (data.e1pt>20.9 && abs(data.e1eta)<2.4) || (data.e2pt>20.9 && abs(data.e2eta)<2.4) ) {
      if (data.mht>87) {
	result+=Trigger::kEleMet;
      }
    }

    if ( (data.m1pt>13.2 && abs(data.m1eta)<2.4) || (data.m2pt>13.2 && abs(data.m2eta)<2.4) ) {
      if (data.mht>57) {
	result+=Trigger::kMuMet;
      }
    }

    if (data.ht>317)
      result+=Trigger::kHt;
    
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
