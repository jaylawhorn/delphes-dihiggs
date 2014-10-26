#ifndef EventShape_h
#define EventShape_h
  
#include "TMatrixDSym.h"
#include "TVectorD.h"
#include "TMath.h"

#include "Math/Vector3D.h"
 
#include <vector>
 
class EventShape {
 
 public:
  EventShape();
  EventShape(const std::vector<ROOT::Math::XYZVector>& inputVectors);  
  ~EventShape(){};

  void addVector(const ROOT::Math::XYZVector& newVector);   

  double isotropy(const unsigned int& numberOfSteps = 1000) const;
   
  double circularity(const unsigned int& numberOfSteps = 1000) const;
  
  double sphericity(double = 2.)  const;

  double aplanarity(double = 2.)  const;

  double C(double = 2.) const;

  double D(double = 2.) const;

  double H(int order) const;

  double ZerothMoment() const;

  double R(int order) const;
   
 private:
  TMatrixDSym compMomentumTensor(double = 2.) const;
  TVectorD compEigenValues(double = 2.) const; 
  std::vector<ROOT::Math::XYZVector> _inputVectors;
};
 
#endif
 
