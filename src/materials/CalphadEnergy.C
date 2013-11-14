/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  1 November 2013
*
*************************************************************************/

#include "CalphadEnergy.h"

template<>
InputParameters validParams<CalphadEnergy>()
{
  InputParameters params = validParams<Material>();

  params.addRequiredParam<std::vector<Real> >("pure_endpoint_low_coeffs", "A, B, C, D, E Gibbs coeffs for low-concentration pure endpoint");
  params.addRequiredParam<std::vector<Real> >("pure_endpoint_high_coeffs", "A, B, C, D, E Gibbs coeffs for high-concentration pure endpoint");
  params.addRequiredParam<std::vector<Real> >("mixture_coeffs", "A, B, C, D, E Gibbs coeffs for mixture endpoint");

  params.addRequiredParam<std::vector<Real> >("L0_coeffs", "Rudlich-Kister L0 polynomial coefficients");
  params.addRequiredParam<std::vector<Real> >("L1_coeffs", "Rudlich-Kister L1 polynomial coefficients");

  params.addRequiredParam<std::vector<Real> >("low_coeffs", "tweak polynomial coeffients for low concentration");
  params.addRequiredParam<std::vector<Real> >("high_coeffs", "tweak polynomial coeffients for high concentration");

  //Default in J/mol-K
  params.addParam<Real>("gas_constant", 8.3144621, "Universal gas constant");

  params.addRequiredCoupledVar("coupled_temperature", "temperature to be used to calculating Gibbs energies");
  params.addRequiredCoupledVar("coupled_concentration", "concentration to be used to calculating Gibbs energies");

  return params;
}

CalphadEnergy::CalphadEnergy(const std::string & name, InputParameters parameters)
    : Material(name, parameters),
      _pure_endpoint_low_coeffs(getParam<std::vector<Real> >("pure_endpoint_low_coeffs")),
      _pure_endpoint_high_coeffs(getParam<std::vector<Real> >("pure_endpoint_high_coeffs")),
      _mixture_coeffs(getParam<std::vector<Real> >("mixture_coeffs")),
      _L0_coeffs(getParam<std::vector<Real> >("L0_coeffs")),
      _L1_coeffs(getParam<std::vector<Real> >("L1_coeffs")),
      _low_coeffs(getParam<std::vector<Real> >("low_coeffs")),
      _high_coeffs(getParam<std::vector<Real> >("high_coeffs")),

      _R(getParam<Real>("gas_constant")),

      _T(coupledValue("coupled_temperature")),
      _c(coupledValue("coupled_concentration"))
{
  //Make sure everything is set up from input correctly
  if( (_pure_endpoint_low_coeffs.size() != 5) ||
      (_pure_endpoint_high_coeffs.size() != 5) ||
      (_mixture_coeffs.size() != 3) ||
      (_L0_coeffs.size() != 2) ||
      (_L1_coeffs.size() != 2) ||
      (_low_coeffs.size() != 10) ||
      (_high_coeffs.size() != 10) )
    mooseError("Please supply the correct # of values for the Gibbs coefficients (CalphadEnergy).");
}

void
CalphadEnergy::computeQpProperties()
{
}

Real
CalphadEnergy::calculateFirstLatticeGminusHser()
{
  return _pure_endpoint_low_coeffs[0]
    + _pure_endpoint_low_coeffs[1]*_T[_qp]
    + _pure_endpoint_low_coeffs[2]*_T[_qp]*std::log(_T[_qp])
    + _pure_endpoint_low_coeffs[3]*_T[_qp]*_T[_qp]
    + _pure_endpoint_low_coeffs[4]/_T[_qp];
}

Real
CalphadEnergy::computeGMix()
{
  return calculateReference() + calculateIdeal() + calculateExcess();
}

Real
CalphadEnergy::computeEndPolynomial(bool low, PolynomialOrder deriv)
{
  if(low)
  {
    switch (deriv)
    {
    case Zero: //calculate low polynomial
      return _low_coeffs[0] + _low_coeffs[1]*_c[_qp] + _low_coeffs[2]*_T[_qp]
        + _low_coeffs[3]*_c[_qp]*_c[_qp] + _low_coeffs[4]*_T[_qp]*_T[_qp]
        + _low_coeffs[5]*_c[_qp]*_c[_qp]*_c[_qp] + _low_coeffs[6]*_T[_qp]*_T[_qp]*_T[_qp]
        + _low_coeffs[7]*_c[_qp]*_T[_qp] + _low_coeffs[8]*_c[_qp]*_c[_qp]*_T[_qp]
        + _low_coeffs[9]*_c[_qp]*_T[_qp]*_T[_qp];

    case First: //calculate low polynomial first derivative
      return  _low_coeffs[1] + 2*_low_coeffs[3]*_c[_qp] + 3*_low_coeffs[5]*_c[_qp]*_c[_qp]
        + _low_coeffs[7]*_T[_qp] + 2*_low_coeffs[8]*_c[_qp]*_T[_qp] + _low_coeffs[9]*_T[_qp]*_T[_qp];

    case Second:  //calculate low polynomial 2nd derivative
      return 2*_low_coeffs[3] + 6*_low_coeffs[5]*_c[_qp] + 2*_low_coeffs[8]*_T[_qp];

    case Third: //calculate low polynomial 3rd derivative
      return  6*_low_coeffs[5];
    }
    mooseError("Invalid type passed in");
  }

  else
  {
    switch (deriv)
    {
    case Zero: //calculate high polynomial
      return _high_coeffs[0] + _high_coeffs[1]*_c[_qp] + _high_coeffs[2]*_T[_qp]
        + _high_coeffs[3]*_c[_qp]*_c[_qp] + _high_coeffs[4]*_T[_qp]*_T[_qp]
        + _high_coeffs[5]*_c[_qp]*_c[_qp]*_c[_qp] + _high_coeffs[6]*_T[_qp]*_T[_qp]*_T[_qp]
        + _high_coeffs[7]*_c[_qp]*_T[_qp] + _high_coeffs[8]*_c[_qp]*_c[_qp]*_T[_qp]
        + _high_coeffs[9]*_c[_qp]*_T[_qp]*_T[_qp];

    case First: //calculate high polynomial first derivative
      return  _high_coeffs[1] + 2*_high_coeffs[3]*_c[_qp] + 3*_high_coeffs[5]*_c[_qp]*_c[_qp]
        + _high_coeffs[7]*_T[_qp] + 2*_high_coeffs[8]*_c[_qp]*_T[_qp] + _high_coeffs[9]*_T[_qp]*_T[_qp];

    case Second:  //calculate high polynomial 2nd derivative
      return 2*_high_coeffs[3] + 6*_high_coeffs[5]*_c[_qp] + 2*_high_coeffs[8]*_T[_qp];

    case Third: //calculate high polynomial 3rd derivative
      return 6*_high_coeffs[5];
    }
       mooseError("Invalid type passed in");
  }
}


Real
CalphadEnergy::calculateReference()
{
  return 0;
}

Real
CalphadEnergy::calculateIdeal()
{
  return 0;
}

Real
CalphadEnergy::calculateExcess()
{
  return 0;
}

Real
CalphadEnergy::calculateSecondLatticeGminusHser()
{
  return 0;
}

Real
CalphadEnergy::computeDGMixDc()
{
  return 0;
}

Real
CalphadEnergy::computeD2GMixDc2()
{
  return 0;
}

Real
CalphadEnergy::computeD3GMixDc3()
{
  return 0;
}
