#include "EventShape.h"
#include "FoxWolframMoments.h"
#include "TMath.h"

using namespace ROOT::Math;  

ClassImp(EventShape)

EventShape::EventShape() 
 {}

EventShape::EventShape(const std::vector<XYZVector>& inputVectors) 
  : _inputVectors(inputVectors)
 {}
 
void 
EventShape::addVector(const XYZVector& newVector)
{
  _inputVectors.push_back(newVector);
}
    
double 
EventShape::isotropy(const unsigned int& numberOfSteps) const
{
  const double deltaPhi=2*TMath::Pi()/numberOfSteps;
  double phi = 0, eIn =-1., eOut=-1.;
  for(unsigned int i=0; i<numberOfSteps; ++i){
    phi+=deltaPhi;
     double sum=0;
     for(unsigned int j=0; j<_inputVectors.size(); ++j){
       // sum over inner product of unit vectors and momenta
       sum+=TMath::Abs(TMath::Cos(phi)*_inputVectors[j].x()+TMath::Sin(phi)*_inputVectors[j].y());
     }
     if( eOut<0. || sum<eOut ) eOut=sum;
     if( eIn <0. || sum>eIn  ) eIn =sum;
  }
  return (eIn-eOut)/eIn;
}
 
double 
EventShape::circularity(const unsigned int& numberOfSteps) const
 {
   const double deltaPhi=2*TMath::Pi()/numberOfSteps;
   double circularity_=-1, phi=0, area = 0;
   for(unsigned int i=0;i<_inputVectors.size();i++) {
     area+=TMath::Sqrt(_inputVectors[i].x()*_inputVectors[i].x()+_inputVectors[i].y()*_inputVectors[i].y());
   }
   for(unsigned int i=0; i<numberOfSteps; ++i){
     phi+=deltaPhi;
     double sum=0, tmp=0.;
     for(unsigned int j=0; j<_inputVectors.size(); ++j){
       sum+=TMath::Abs(TMath::Cos(phi)*_inputVectors[j].x()+TMath::Sin(phi)*_inputVectors[j].y());
     }
     tmp=TMath::Pi()/2*sum/area;
     if( circularity_<0 || tmp<circularity_ ){
       circularity_=tmp;
     }
   }
   return circularity_;
 }
 
TMatrixDSym 
EventShape::compMomentumTensor(double r) const
 {
   TMatrixDSym momentumTensor(3);
   momentumTensor.Zero();
 
   if ( _inputVectors.size() < 2 ){
     return momentumTensor;
   }
 
   // fill momentumTensor from inputVectors
   double norm = 1.;
   for ( int i = 0; i < (int)_inputVectors.size(); ++i ){
     double p2 = _inputVectors[i].Dot(_inputVectors[i]);
     double pR = ( r == 2. ) ? p2 : TMath::Power(p2, 0.5*r);
     norm += pR;
     double pRminus2 = ( r == 2. ) ? 1. : TMath::Power(p2, 0.5*r - 1.);
     momentumTensor(0,0) += pRminus2*_inputVectors[i].x()*_inputVectors[i].x();
     momentumTensor(0,1) += pRminus2*_inputVectors[i].x()*_inputVectors[i].y();
     momentumTensor(0,2) += pRminus2*_inputVectors[i].x()*_inputVectors[i].z();
     momentumTensor(1,0) += pRminus2*_inputVectors[i].y()*_inputVectors[i].x();
     momentumTensor(1,1) += pRminus2*_inputVectors[i].y()*_inputVectors[i].y();
     momentumTensor(1,2) += pRminus2*_inputVectors[i].y()*_inputVectors[i].z();
     momentumTensor(2,0) += pRminus2*_inputVectors[i].z()*_inputVectors[i].x();
     momentumTensor(2,1) += pRminus2*_inputVectors[i].z()*_inputVectors[i].y();
     momentumTensor(2,2) += pRminus2*_inputVectors[i].z()*_inputVectors[i].z();
   }
 
   //std::cout << "momentumTensor:" << std::endl;
   //std::cout << momentumTensor(0,0) << " " << momentumTensor(0,1) << " " << momentumTensor(0,2) 
   //          << momentumTensor(1,0) << " " << momentumTensor(1,1) << " " << momentumTensor(1,2) 
   //          << momentumTensor(2,0) << " " << momentumTensor(2,1) << " " << momentumTensor(2,2) << std::endl;
 
   // return momentumTensor normalized to determinant 1
   return (1./norm)*momentumTensor;
 }
 
TVectorD
EventShape::compEigenValues(double r) const
{
  TVectorD eigenValues(3);
  TMatrixDSym myTensor = compMomentumTensor(r);
  if( myTensor.IsSymmetric() ){
    if( myTensor.NonZeros() != 0 ) myTensor.EigenVectors(eigenValues);
  }
  
   // CV: TMatrixDSym::EigenVectors returns eigen-values and eigen-vectors
   //     ordered by descending eigen-values, so no need to do any sorting here...
   //std::cout << "eigenValues(0) = " << eigenValues(0) << ","
   //          << " eigenValues(1) = " << eigenValues(1) << ","
   //          << " eigenValues(2) = " << eigenValues(2) << std::endl;
 
  return eigenValues;
}
 
double 
EventShape::sphericity(double r) const
{
  TVectorD eigenValues = compEigenValues(r);
  return 1.5*(eigenValues(1) + eigenValues(2));
}

double 
 EventShape::aplanarity(double r) const
{
  TVectorD eigenValues = compEigenValues(r);
  return 1.5*eigenValues(2);
}
 
double 
EventShape::C(double r) const
{
  TVectorD eigenValues = compEigenValues(r);
   return 3.*(eigenValues(0)*eigenValues(1) + eigenValues(0)*eigenValues(2) + eigenValues(1)*eigenValues(2));
}
 
double 
EventShape::D(double r) const
{
  TVectorD eigenValues = compEigenValues(r);
  return 27.*eigenValues(0)*eigenValues(1)*eigenValues(2);
}
 
double 
EventShape::H(int i) const
{
  FoxWolframMoments fwm(i);
  fwm.Compute(_inputVectors);
  return fwm.H(i);
}

double 
EventShape::ZerothMoment() const
{
  FoxWolframMoments fwm(0);
  fwm.Compute(_inputVectors);
  return fwm.ZerothMoment();
}
 
double 
EventShape::R(int i) const
{
  FoxWolframMoments fwm(i);
  fwm.Compute(_inputVectors);
  return fwm.R(i);
}

 
