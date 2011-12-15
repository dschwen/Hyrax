/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  16 November 2011
*
*  This code inherits from CHBulk in ELK
*  
*  This code handles the materials parameters for a coupled 
*  conserved order parameter, non-conserved order parameter
*  system. 
*************************************************************************/


#ifndef PFMobilityLandau_H
#define PFMobilityLandau_H

#include "Material.h"

//Forward Declarations
class PFMobilityLandau;

template<>
InputParameters validParams<PFMobilityLandau>();

class PFMobilityLandau : public Material
{
public:
  PFMobilityLandau(const std::string & name,
          InputParameters parameters);
  
protected:
  virtual void computeProperties();

private:
  
  // Variables for getting values from input file
  Real _mob_CH;
  Real _kappa_CH;

  Real _mob_AC;
  Real _kappa_AC;  

  Real _a1_i;
  Real _a2_i;
  Real _a3_i;
  Real _a4_i;

  Real _c1_i;
  Real _c2_i;

  // Cahn-Hilliard equation 
  MaterialProperty<Real> & _M;
  MaterialProperty<RealGradient> & _grad_M;
  MaterialProperty<Real> & _kappa_c;

  // Alan-Cahn equation 
  MaterialProperty<Real> & _L;
  MaterialProperty<Real> & _kappa_n;

  // Landau polynomial coefficients
  MaterialProperty<Real> & _a1;
  MaterialProperty<Real> & _a2;
  MaterialProperty<Real> & _a3;
  MaterialProperty<Real> & _a4;

  // Landau polynomial parameters
  MaterialProperty<Real> & _c1;
  MaterialProperty<Real> & _c2;

};

#endif //PFMobilityLandau_H
