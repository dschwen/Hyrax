#include "CalphadAB1CD2Old.h"
//#include "CalphadParameters.h"

#include <cmath>

CalphadAB1CD2Old::CalphadAB1CD2Old() :
    CalphadFreeEnergy(),
    _pure_EP1_phase1_coeffs()
{
}


void
CalphadAB1CD2Old::parameterize(Real & R, std::vector<Real> & low, std::vector<Real> & high, std::vector<Real> & mix,
                            std::vector<Real> L0, std::vector<Real> L1, std::vector<Real> Pure1)
{
  _pure_endpoint_1_coeffs = low;
  _pure_endpoint_2_coeffs = high;
  _mixture_coeffs = mix;

  _L0_coeffs = L0;
  _L1_coeffs = L1;

  _pure_EP1_phase1_coeffs = Pure1;

  _R = R;
}

Real
CalphadAB1CD2Old::calculateReference(const Real & c, const Real & T) const
{
  Real first_energy = calculateFirstLatticeGminusHser(c, T);
  Real second_energy = calculateSecondLatticeGminusHser(c, T);

  return  0.5*( (2 - 3*c)*first_energy + c*second_energy );
}


Real
CalphadAB1CD2Old::calculateIdeal(const Real & c, const Real & T) const
{
  return _R*T*( c*std::log(c/(2-2*c)) + (2-3*c)*std::log((2-3*c)/(2-2*c)) );
}

Real
CalphadAB1CD2Old::calculateExcess(const Real & c, const Real & T) const
{
  Real L0 = _L0_coeffs[0] + _L0_coeffs[1]*T;
  Real L1 = _L1_coeffs[0] + _L1_coeffs[1]*T;

  Real polynomial = (c-1)*L0 + (1-2*c)*L1;

  Real prefactor = 3*c*c - 2*c;
  prefactor /= 4*(c-1)*(c-1);

  return prefactor*polynomial;
}

Real
CalphadAB1CD2Old::calculateSecondLatticeGminusHser(const Real & c, const Real & T) const
{
  Real first_term = _mixture_coeffs[0]
    + _mixture_coeffs[1]*T
    + _mixture_coeffs[2]*T*std::log(T);

  Real second_term = _pure_EP1_phase1_coeffs[0]
    + _pure_EP1_phase1_coeffs[1]*T
    + _pure_EP1_phase1_coeffs[2]*T*std::log(T)
    + _pure_EP1_phase1_coeffs[3]*T*T
    + _pure_EP1_phase1_coeffs[4]/T;

  Real third_term = _pure_endpoint_2_coeffs[0]
    + _pure_endpoint_2_coeffs[1]*T
    + _pure_endpoint_2_coeffs[2]*T*std::log(T)
    + _pure_endpoint_2_coeffs[3]*T*T
    + _pure_endpoint_2_coeffs[4]/T;

  return first_term + second_term + third_term;
}

Real
CalphadAB1CD2Old::computeGMix(const Real & c, const Real & T) const
{
  Real c1;
  //make this piecewise in concentration space
  if( c < 0.001)
  {
    c1 = 0.001;

    return CalphadFreeEnergy::computeGMix(c1, T) + computeDGMixDc(c1, T)*(c - c1);
  }
  else if (c > 0.655)
  {
    c1 = 0.655;

    return CalphadFreeEnergy::computeGMix(c1, T) + computeDGMixDc(c1, T)*(c - c1);
  }
  else
    return CalphadFreeEnergy::computeGMix(c, T);
}


Real
CalphadAB1CD2Old::computeDGMixDc(const Real & c, const Real & T) const
{
  Real ref;
  Real ideal;
  Real c1;

 if( c < 0.001)
   c1 = 0.001;

 else if (c > 0.655)
   c1 = 0.655;

 else
   c1 = c;

 ref = -1.5*calculateFirstLatticeGminusHser(c1, T) + 0.5*calculateSecondLatticeGminusHser(c1, T);

  ideal = _R*T*( std::log(4*c1/(1-c1)) - 3*std::log((2-3*c1)/(1-c1)) );

  Real xs;

  Real L0 = _L0_coeffs[0] + _L0_coeffs[1]*T;
  Real L1 = _L1_coeffs[0] + _L1_coeffs[1]*T;

  Real polynomial = (3*c1*c1*c1 -9*c1*c1 + 8*c1 -2)*L0 + (-6*c1*c1*c1 + 18*c1*c1 -12*c1 +2)*L1;

  Real prefactor = 1/(4*(c1-1)*(c1-1)*(c1-1) );

  xs = prefactor*polynomial;

  return ref + ideal + xs;


}

Real
CalphadAB1CD2Old::computeD2GMixDc2(const Real & c, const Real & T) const
{
  Real ref;
  Real ideal;
  Real xs;

  //make this piecewise in concentration space
  if( c < 0.001)
    return 0;

  else if (c > 0.655)
    return 0;

  else
  {
    ref = 0;

    ideal = _R*T*( 2/( c*(3*c*c - 5*c + 2)) );

    Real L0 = _L0_coeffs[0] + _L0_coeffs[1]*T;
    Real L1 = _L1_coeffs[0] + _L1_coeffs[1]*T;

    Real polynomial = (c - 1)*L0 + (3 - 6*c)*L1;

    Real prefactor = 0.5/(4*std::pow(c - 1, 4));

    xs = prefactor*polynomial;

    return ref + ideal + xs;
  }
}
