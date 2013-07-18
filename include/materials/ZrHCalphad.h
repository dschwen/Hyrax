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
  virtual void computeQpProperties();



  Real computeGhcpZr();
  Real computeGhcpZrH();
  Real computeGfccZr();
  Real computeGfccZrH2();
  Real computeGH2();

  Real computeMobility();
  RealGradient computeGradMobility();

private:

  //INPUT VARIABLES
  //Molar Gibbs free energy coefficents
  std::vector<Real> _V_hcp_Zr; //"V" for "victory"...I mean, "vector"!
  std::vector<Real> _V_fcc_Zr;
  std::vector<Real> _V_hcp_ZrH;
  std::vector<Real> _V_fcc_ZrH2;
  std::vector<Real> _V_H2;
  std::vector<Real> _V_L0;
  std::vector<Real> _V_L1;

  //Diffusion coefficient information
  Real _H_Zr_D0;
  Real _H_ZrH2_D0;
  Real _H_Zr_Q0;
  Real _H_ZrH2_Q0;

  Real _mobility_AC;
  Real _kappa_CH;
  Real _kappa_AC;

  Real _T;
  Real _R;

  unsigned int _n_OP_variables;

  ///MATERIALS PROPERTIES
  //molar Gibbs free energies
  MaterialProperty<Real> & _G_hcp_Zr;
  MaterialProperty<Real> & _G_fcc_Zr;
  MaterialProperty<Real> & _G_hcp_ZrH;
  MaterialProperty<Real> & _G_fcc_ZrH2;
  MaterialProperty<Real> & _G_H2;
  MaterialProperty<Real> & _L0;
  MaterialProperty<Real> & _L1;

  MaterialProperty<Real> & _Zr_molar_volume;
  MaterialProperty<Real> & _ZrH2_molar_volume;

  MaterialProperty<Real> & _M;			//Cahn-Hilliard mobility (isotropic)
  MaterialProperty<RealGradient> & _grad_M;
  MaterialProperty<Real> & _kappa_c;		//CH gradient energy coefficient (isotropic)
  MaterialProperty<Real> & _L;			//Allen-Cahn kinetic coefficient (isotropic)
  MaterialProperty<Real> & _kappa_n;		//AC gradient energy coefficient (isotropic)

  // INTERNAL VARIABLES
  Real _D_alpha;                                //Diffusion coefficient of H in alpha Zr
  Real _D_delta;                                //Diffusion coefficient of H in delta ZrH2;

  VariableValue & _X; // conserved variable (atomic fraction)
  VariableGradient & _grad_X;
  std::vector<VariableValue *> _OP; //nonconserved variables (structural order parameter)
  std::vector<VariableGradient *> _grad_OP;



};

#endif //ZRHCALPHAD_H
