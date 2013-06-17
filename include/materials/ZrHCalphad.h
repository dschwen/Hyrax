/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  13 June 2013
*
*************************************************************************/

#ifndef ZRHCALPHAD_H
#define ZRHCALPHAD_H

#include "Material.h"

//forward declaration
class ZrHCalphad;

template<>
InputParameters validParams<ZrHCalphad>();

class ZrHCalphad : public Material
{
public:
  ZrHCalphad(const std::string & name, InputParameters parameters);

protected:
  virtual void computeProperties();

private:

  //INPUT VARIABLES
  std::vector<Real> _G_hcp_Zr_vector;
  std::vector<Real> _G_fcc_Zr_vector;
  std::vector<Real> _G_hcp_ZrH_vector;
  std::vector<Real> _G_fcc_ZrH2_vector;
  std::vector<Real> _G_H2_vector;
  std::vector<Real> _L0_vector;
  std::vector<Real> _L1_vector;

  Real _H_coeff;
  Real _kappa_CH;

  Real _mobility_AC;
  Real _kappa_AC;

  //MATERIALS PROPERTIES
  //molar Gibbs free energy coefficients for HCP Zr
  MaterialProperty<Real> & _G_hcp_Zr_a;
  MaterialProperty<Real> & _G_hcp_Zr_b;
  MaterialProperty<Real> & _G_hcp_Zr_c;
  MaterialProperty<Real> & _G_hcp_Zr_d;
  MaterialProperty<Real> & _G_hcp_Zr_e;

  //molar Gibbs free energy coefficients for FCC Zr
  MaterialProperty<Real> & _G_fcc_Zr_a;
  MaterialProperty<Real> & _G_fcc_Zr_b;
  MaterialProperty<Real> & _G_fcc_Zr_c;
  MaterialProperty<Real> & _G_fcc_Zr_d;
  MaterialProperty<Real> & _G_fcc_Zr_e;

  //molar Gibbs free energy coefficients for HCP ZrH
  MaterialProperty<Real> & _G_hcp_ZrH_a;
  MaterialProperty<Real> & _G_hcp_ZrH_b;

  //molar Gibbs free energy coefficients for FCC ZrH2
  MaterialProperty<Real> & _G_fcc_ZrH2_a;
  MaterialProperty<Real> & _G_fcc_ZrH2_b;
  MaterialProperty<Real> & _G_fcc_ZrH2_c;

  //molar Gibbs free energy coefficients for H2 gas
  MaterialProperty<Real> & _G_H2_a;
  MaterialProperty<Real> & _G_H2_b;
  MaterialProperty<Real> & _G_H2_c;
  MaterialProperty<Real> & _G_H2_d;
  MaterialProperty<Real> & _G_H2_e;

  //Redlich-Kister polynomial coefficients
  MaterialProperty<Real> & _L0_a;
  MaterialProperty<Real> & _L0_b;
  MaterialProperty<Real> & _L1_a;
  MaterialProperty<Real> & _L1_b;


  MaterialProperty<Real> & _M;			//Cahn-Hilliard mobility (isotropic)
  MaterialProperty<RealGradient> & _grad_M;
  MaterialProperty<Real> & _kappa_c;		//CH gradient energy coefficient (isotropic)

  MaterialProperty<Real> & _L;			//Allen-Cahn kinetic coefficient (isotropic)
  MaterialProperty<Real> & _kappa_n;		//AC gradient energy coefficient (isotropic)

};

#endif //ZRHCALPHAD_H
