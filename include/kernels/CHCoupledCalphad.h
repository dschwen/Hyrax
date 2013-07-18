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
  void computeDHeaviside();

  RealGradient computeGradConservedTerm();
  RealGradient computeGradNonconservedTerm();

  Real computeDGalphaDx();
  Real computeDGdeltaDx();
  Real computeD2GalphaDx2();
  Real computeD2GdeltaDx2();
  Real computeD3GalphaDx3();
  Real computeD3GdeltaDx3();

private:
  //molar Gibbs free energies
  MaterialProperty<Real> & _G_hcp_Zr;
  MaterialProperty<Real> & _G_hcp_ZrH;
  MaterialProperty<Real> & _G_fcc_Zr;
  MaterialProperty<Real> & _G_fcc_ZrH2;
  MaterialProperty<Real> & _G_H2;
  MaterialProperty<Real> & _L0;
  MaterialProperty<Real> & _L1;

  MaterialProperty<Real> & _molarVol_alpha_Zr;
  MaterialProperty<Real> & _molarVol_delta_ZrH2;


  Real _T; //temperature
  Real _R; //universal gas constant

  unsigned int _n_OP_variables;
  std::vector<VariableValue *> _OP;
  std::vector<VariableGradient *> _grad_OP;

  Real _Heaviside;
  std::vector<Real> _dHeaviside;
};

#endif //CHCOUPLEDCALPHAD_H
