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
//  std::cout<<"in CalphadAB1CD1 computeQpProperties"<<std::endl;
  _G_AB1CD1[_qp] = computeGMix();
  _dG_dc[_qp] = computeDGMixDc();
  _d2G_dc2[_qp] = computeD2GMixDc2();
  _d3G_dc3[_qp] = computeD3GMixDc3();

//  std::cout<<"end CalphadAB1CD1 computeQpProperties"<<std::endl;
}

Real
CalphadAB1CD1::calculateReference()
{
  Real first_energy = calculateFirstLatticeGminusHser();
  Real second_energy = calculateSecondLatticeGminusHser();

  return  (1 - 2*_c[_qp])*first_energy + _c[_qp]*second_energy;
}


Real
CalphadAB1CD1::calculateIdeal()
{

 return _R*_T[_qp]*( _c[_qp]*std::log(_c[_qp]/(1-_c[_qp]))
                      + (1-2*_c[_qp])*std::log((1-2*_c[_qp])/(1-_c[_qp])) );
}

Real
CalphadAB1CD1::calculateExcess()
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
CalphadAB1CD1::computeGMix()
{
  //make this piecewise in concentration space
  if( _c[_qp] < 0.0017)
  {
    return computeEndPolynomial(true, Zero);
  }
  else if (_c[_qp] > 0.49)
  {
    return  computeEndPolynomial(false, Zero);

  }

  else
  {

    return CalphadEnergy::computeGMix();
  }
}

Real
CalphadAB1CD1::computeDGMixDc()
{
  Real ref;
  Real ideal;

  //make this piecewise in concentration space
  if( _c[_qp] < 0.0017)
  {
    return computeEndPolynomial(true, First);
  }
  else if (_c[_qp] > 0.49)
  {
    return computeEndPolynomial(false, First);
  }

  else
  {
    ref = -2*calculateFirstLatticeGminusHser() + calculateSecondLatticeGminusHser();

    ideal = _R*_T[_qp]*( std::log(_c[_qp]/(1-_c[_qp])) - 2*std::log((1-2*_c[_qp])/(1-_c[_qp])) );

    return ref + ideal;
  }

}

Real
CalphadAB1CD1::computeD2GMixDc2()
{
  Real ref;
  Real ideal;

  //make this piecewise in concentration space
  if( _c[_qp] < 0.0017)
  {
    return computeEndPolynomial(true, Second);
  }
  else if (_c[_qp] > 0.49)
  {
    return computeEndPolynomial(false, Second);
  }

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
  if( _c[_qp] < 0.0017)
  {
    return computeEndPolynomial(true, Third);
  }
  else if (_c[_qp] > 0.49)
  {
    return computeEndPolynomial(false, Third);
  }

  else
  {
    ref = 0;

    ideal = -1*_R*_T[_qp]*( (6*_c[_qp]*_c[_qp] -6*_c[_qp] + 1)
                            / std::pow( _c[_qp]*(2*_c[_qp]*_c[_qp] - 3*_c[_qp] + 1), 2) );

    return ref + ideal;
  }

}

