/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  7 November 2013
*
*************************************************************************/

#include "CalphadAB1CD1.h"

template<>
InputParameters validParams<CalphadAB1CD1>()
{
  InputParameters params = validParams<CalphadEnergy>();

  return params;
}

CalphadAB1CD1::CalphadAB1CD1(const std::string & name, InputParameters parameters) :
    CalphadEnergy(name, parameters),
    _G_AB1CD1(declareProperty<Real>("G_AB1CD1")),
    _dG_dc(declareProperty<Real>("dGAB1CD1_dc")),
    _d2G_dc2(declareProperty<Real>("d2GAB1CD1_dc2")),
    _d3G_dc3(declareProperty<Real>("d3GAB1CD1_dc3"))
{
}

void
CalphadAB1CD1::computeQpProperties()
{
  Real c;

  if (_c[_qp]  < 0.001)
    c = 0.001;
  else if (_c[_qp] > 0.499)
    c = 0.499;
  else c = _c[_qp];

  _G_AB1CD1[_qp] = computeGMix(c);
  _dG_dc[_qp] = computeDGMixDc(c);
  _d2G_dc2[_qp] = computeD2GMixDc2();
  _d3G_dc3[_qp] = computeD3GMixDc3();

}

Real
CalphadAB1CD1::calculateReference(Real c)
{
  Real first_energy = calculateFirstLatticeGminusHser();
  Real second_energy = calculateSecondLatticeGminusHser();

  return  (1 - 2*c)*first_energy + c*second_energy;
}


Real
CalphadAB1CD1::calculateIdeal(Real c)
{
  return _R*_T[_qp]*( c*std::log(c/(1-c)) + (1-2*c)*std::log((1-2*c)/(1-c)) );
}

Real
CalphadAB1CD1::calculateExcess(Real c)
{
  return 0;
}

Real
CalphadAB1CD1::calculateSecondLatticeGminusHser()
{
  Real first_term = _mixture_coeffs[0] + _mixture_coeffs[1]*_T[_qp];

  Real second_term = _pure_endpoint_low_coeffs[0]
    + _pure_endpoint_low_coeffs[1]*_T[_qp]
    + _pure_endpoint_low_coeffs[2]*_T[_qp]*std::log(_T[_qp])
    + _pure_endpoint_low_coeffs[3]*_T[_qp]*_T[_qp]
    + _pure_endpoint_low_coeffs[4]/_T[_qp];

  Real third_term = 0.5*(_pure_endpoint_high_coeffs[0]
                           + _pure_endpoint_high_coeffs[1]*_T[_qp]
                           + _pure_endpoint_high_coeffs[2]*_T[_qp]*std::log(_T[_qp])
                           + _pure_endpoint_high_coeffs[3]*_T[_qp]*_T[_qp]
                           + _pure_endpoint_high_coeffs[4]/_T[_qp] );

  return first_term + second_term + third_term;
}

Real
CalphadAB1CD1::computeGMix(Real c)
{
  Real c1;
  //make this piecewise in concentration space
 if(_c[_qp] < 0.001)
  {
    c1 = 0.001;

    return CalphadEnergy::computeGMix(c1) + computeDGMixDc(c1)*(_c[_qp] - c1);
  }

  else if (_c[_qp] > 0.499)
  {
    c1 = 0.499;

    return CalphadEnergy::computeGMix(c1) + computeDGMixDc(c1)*(_c[_qp] - c1);
  }

  else
    return CalphadEnergy::computeGMix(_c[_qp]);
}

Real
CalphadAB1CD1::computeDGMixDc(Real c)
{
  Real ref;
  Real ideal;

  ref = -2*calculateFirstLatticeGminusHser() + calculateSecondLatticeGminusHser();
  ideal = _R*_T[_qp]*( std::log(c/(1-c)) - 2*std::log((1-2*c)/(1-c)) );

  return ref + ideal;
}

Real
CalphadAB1CD1::computeD2GMixDc2()
{
  Real ref;
  Real ideal;

  //make this piecewise in concentration space
  if( _c[_qp] < 0.001)
    return 0;

  else if (_c[_qp] > 0.499)
    return 0;

  else
  {
    ref = 0;

    ideal = _R*_T[_qp]*( 1/( _c[_qp]*(2*_c[_qp]*_c[_qp] - 3*_c[_qp] + 1)) );

    return ref + ideal;
  }
}

Real
CalphadAB1CD1::computeD3GMixDc3()
{
  Real ref;
  Real ideal;

  //make this piecewise in concentration space
  if( _c[_qp] < 0.001)
    return 0;

  else if (_c[_qp] > 0.499)
    return 0;

  else
  {
    ref = 0;

    ideal = -1*_R*_T[_qp]*( (6*_c[_qp]*_c[_qp] -6*_c[_qp] + 1)
                            / std::pow( _c[_qp]*(2*_c[_qp]*_c[_qp] - 3*_c[_qp] + 1), 2) );

    return ref + ideal;
  }
}

