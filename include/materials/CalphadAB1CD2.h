/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  7 November 2013
*
*************************************************************************/

#ifndef CALPHADAB1CD2_H
#define CALPHADAB1CD2_H

#include "CalphadEnergy.h"

class CalphadAB1CD2;

template<>
InputParameters validParams<CalphadAB1CD2>();

class CalphadAB1CD2 : public CalphadEnergy
{
public:
  CalphadAB1CD2(const std::string & name, InputParameters parameters);

protected:
  virtual void computeQpProperties();

  virtual Real calculateReference(Real c);
  virtual Real calculateIdeal(Real c);
  virtual Real calculateExcess(Real c);

  virtual Real calculateSecondLatticeGminusHser();

  virtual Real computeGMix(Real c);
  virtual Real computeDGMixDc(Real c);
  virtual Real computeD2GMixDc2();
  virtual Real computeD3GMixDc3();

private:

  std::vector<double> _pure_EP1_phase1_coeffs;
  MaterialProperty<Real> & _G_AB1CD2;
  MaterialProperty<Real> & _dG_dc;
  MaterialProperty<Real> & _d2G_dc2;
  MaterialProperty<Real> & _d3G_dc3;
};

#endif //CALPHADAB1CD2_H
