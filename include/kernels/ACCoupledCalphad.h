/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  12 June 2013
*
*************************************************************************/

#ifndef ACCOUPLEDCALPHAD_H
#define ACCOUPLEDCALPHAD_H

#include "ACBulk.h"

//Forward declarations
class ACCoupledCalphad;

template<>
InputParameters validParams<ACCoupledCalphad>();

class ACCoupledCalphad : public ACBulk
{
public:
  ACCoupledCalphad(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeDFDOP(PFFunctionType type);

  Real computeGalphamix();
  Real computeGdeltamix();

  Real computeGhcpZr();
  Real computeGhcpZrH();
  Real computeGfccZr();
  Real computeGfccZrH2();
  Real computeGH2();

private:
  //molar Gibbs free energy coefficients for HCP Zr
  MaterialProperty<Real> & _G_hcp_Zr_a;
  MaterialProperty<Real> & _G_hcp_Zr_b;
  MaterialProperty<Real> & _G_hcp_Zr_c;
  MaterialProperty<Real> & _G_hcp_Zr_d;
  MaterialProperty<Real> & _G_hcp_Zr_e;

  //molar Gibbs free energy coefficients for HCP ZrH
  MaterialProperty<Real> & _G_hcp_ZrH_a;
  MaterialProperty<Real> & _G_hcp_ZrH_b;

  //molar Gibbs free energy coefficients for FCC Zr
  MaterialProperty<Real> & _G_fcc_Zr_a;
  MaterialProperty<Real> & _G_fcc_Zr_b;
  MaterialProperty<Real> & _G_fcc_Zr_c;
  MaterialProperty<Real> & _G_fcc_Zr_d;
  MaterialProperty<Real> & _G_fcc_Zr_e;

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

  Real _T; //temperature
  Real _R; //universal gas constant

  Real _w; //free energy well height

  VariableValue & _coupled_CH_var;  //Cahn-Hilliard (conserved) variable

  unsigned int _n_OP_vars;
  unsigned int _OP_number;

  std::vector<VariableValue *> _coupled_OP_vars;
};

#endif //ACCOUPLEDCALPHAD_H
