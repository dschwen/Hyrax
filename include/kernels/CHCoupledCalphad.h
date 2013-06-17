/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  11 June 2013
*
*************************************************************************/

#ifndef CHCOUPLEDCALPHAD_H
#define CHCOUPLEDCALPHAD_H

#include "CHBulk.h"

/**
 * CHCoupledCalphad handles the conserved order parameter (atomic fraction)
 * evolved using the Cahn-Hilliard equation coupled to an arbitrary number of
 * non-conserved structural order parameters for the Zr-delta hydride system.
 * It uses free energy equations taken from A.T. Dinsdale 15 (1991) CALPHAD 317,
 * N. Dupin et al, J. Nuc. Mater. 275 (1999) 287 and J.D. Cox et al., CODATA Key
 * Values for Thermodynamics, Hemisphere Publishing Corporation, 1989.
 */

//Forward declarations
class CHCoupledCalphad;

template<>
InputParameters validParams<CHCoupledCalphad>();

class CHCoupledCalphad : public CHBulk
{
public:
  CHCoupledCalphad(const std::string & name, InputParameters parameters);

protected:
  virtual RealGradient computeGradDFDCons(PFFunctionType type, Real c, RealGradient grad_c);

  Real computeHeaviside();
  RealGradient computeGradConservedTerm(Real & h);
  RealGradient computeGradNonconservedTerm();

  Real computeDGalphaDx();
  Real computeDGdeltaDx();

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

  unsigned int _n_OP_variables;
  std::vector<VariableValue *> _coupled_OP_vars;
  std::vector<VariableGradient *> _coupled_OP_grads;
};

#endif //CHCOUPLEDCALPHAD_H
