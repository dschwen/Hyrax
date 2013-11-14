/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  1 November 2013
*
*************************************************************************/

#ifndef CALPHADENERGY_H
#define CALPHADENERGY_H

#include "Material.h"

//forward declaration
class CalphadEnergy;

template<>
InputParameters validParams<CalphadEnergy>();

class CalphadEnergy : public Material
{
public:
  CalphadEnergy(const std::string & name, InputParameters parameters);

protected:
  enum PolynomialOrder
  {
    Zero,
    First,
    Second,
    Third
  };

  virtual void computeQpProperties();

  virtual Real calculateReference();
  virtual Real calculateIdeal();
  virtual Real calculateExcess();
  virtual Real calculateFirstLatticeGminusHser();
  virtual Real calculateSecondLatticeGminusHser();

  virtual Real computeGMix();
  virtual Real computeDGMixDc();
  virtual Real computeD2GMixDc2();
  virtual Real computeD3GMixDc3();
  virtual Real computeEndPolynomial(bool low, PolynomialOrder deriv);

  //vectors to hold the coefficients of the 2 endpoints for mixing
  std::vector<Real> _pure_endpoint_low_coeffs;
  std::vector<Real> _pure_endpoint_high_coeffs;
  std::vector<Real> _mixture_coeffs;

  //vectors to hold the Rudlich-Kister polynomial coefficents
  std::vector<Real> _L0_coeffs;
  std::vector<Real> _L1_coeffs;

  //vectors to hold the coefficients of the calphad endpoint polynomial tweak
  std::vector<Real> _low_coeffs;
  std::vector<Real> _high_coeffs;

  Real _R;                                      //Universal gas constant

  //COUPLED VARIABLES
  VariableValue & _T;                           //coupled Temperature field
  VariableValue & _c;                           //coupled concentration field

private:

  // Inherited classes will put their energies as materials properties
  // here and declare with unique names

};

#endif //CALPHADENERGY_H
