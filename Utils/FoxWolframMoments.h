#ifndef FWM_H
#define FWM_H

#include "TVector3.h"
#include "TVector.h"

#include "Math/Vector3D.h"

class FoxWolframMoments {

 public:
  FoxWolframMoments( int order=4 );
  ~FoxWolframMoments(){};

  void Compute(const std::vector<ROOT::Math::XYZVector>& inputVectors);
  void Reset();

  const TVector& Moments() const      { return _FWarray; }   
  const TVector& SumArray() const     { return _sumarray; } 
  double H(int order)  const     { return _FWarray(order); }
  double ZerothMoment() const     { return _FWarray(0); } 
  double R( int order ) const;  // normalized to zeroth-moment

  static double Legendre( int l, int m, double x );
 
private:
  int _nmom;
  TVector _FWarray;
  TVector _sumarray;

};

#endif
