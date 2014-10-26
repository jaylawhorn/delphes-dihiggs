#include "FoxWolframMoments.h"
#include <TMath.h>
#include <assert.h>
	
using namespace ROOT::Math;  
	
ClassImp(FoxWolframMoments)
		
FoxWolframMoments::FoxWolframMoments( int maxorder )
  : _nmom( maxorder+1 )
  , _FWarray( _nmom )
  , _sumarray( _nmom )
{}

void
FoxWolframMoments::Compute(const std::vector<XYZVector>& inputVectors)
{
  // initialize
  Reset();
	 
  if( inputVectors.size()==0 ) return;
	
  double s = 0.;
  int l;
  
  // start a loop over the all candidates
  for(unsigned int i=0; i<inputVectors.size(); ++i){
    // this candidate's 3-momentum
    TVector3 p1(inputVectors[i].X(),inputVectors[i].Y(),inputVectors[i].Z());  // ok, this is somewhat ...
    double pmag1 = p1.Mag();
    
    // loop over other candidates, starting at the next one in the list
      for(unsigned int j=i; j<inputVectors.size(); ++j){
	// this candidate's 3-momentum
	TVector3 p2(inputVectors[j].X(),inputVectors[j].Y(),inputVectors[j].Z());
	double pmag2 = p2.Mag();
	  
	// the cosine of the angle between the two candidates
	double cosPhi =  cos ( p1.Angle(p2) );
	
	// the contribution of this pair of track
	// (note the factor 2 : the pair enters the sum twice)
	for( l=0; l<_nmom; l++ )
	  _sumarray(l) += 2 * pmag1 * pmag2 * Legendre( l, 0, cosPhi );
      }
      
      // contribution of this track
      for( l=0; l<_nmom; l++ )
	_sumarray(l) += pmag1 * pmag1 * Legendre( l, 0, 1. );
      
      // total energy
      s += p1.Mag();
      
    }
  
  // well ...
  if( s<=0. ) return;
	 
  // normalize Fox Wolfram Moments
  for(int i=0; i<_nmom; i++)
    _FWarray(i) = _sumarray(i)/pow(s,2) ;
  
}
	
void
FoxWolframMoments::Reset() 
{
  for ( int i=0; i<_nmom; i++) 
    { 
      _FWarray(i) = 0.;
      _sumarray(i) = 0.;
    }
}
	
double
FoxWolframMoments::R( int order ) const
{
  if( H(0)>0. ) {
    if( order < _nmom ){ 
      return ( H(order)/H(0) );
    }
  }
  return 0.;
}

double 
FoxWolframMoments::Legendre( int l, int m, double x )
{
  assert(m >= 0.);
  assert(m <= l);
  assert(fabs(x) <= 1.);
  
  double pmm = 1.;
  
  if(m > 0)
    {
      double somx2 = sqrt((1. - x) * (1. + x));
      double fact = 1.;
      
      for(int i=0; i<m; i++)
	{
	  pmm *= -fact * somx2;
	  fact += 2.0;
	}
    }
	
  if(l == m)
    return pmm;
  
  else
    {
      double pmmp1 = x * (2 * m + 1) * pmm;
      if(l == m + 1)
	return pmmp1;
      else
	{
	  for(int ll=m+2; ll<=l; ll++)
	            {
		      double pll = (x * (2 * ll - 1) * pmmp1 - (ll + m - 1) * pmm) / 
			(ll - m);
		      pmm = pmmp1;
	                pmmp1 = pll;
	            }
	
	  return pmmp1;
	}
    }
}

