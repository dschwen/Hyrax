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

  Real computeDHeavisideDOP();
  Real computeD2HeavisideDOP2();
  Real computeDBarrierDOP(Real & SS, Real & QS, Real & SM);
  Real computeD2BarrierDOP2(Real & SS, Real & QS, Real & SM);

  Real computeGalphamix();
  Real computeGdeltamix();

private:
  //molar Gibbs free energy coefficients
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

  Real _w; //free energy well height

  VariableValue & _X;  //Cahn-Hilliard (conserved) variable

  unsigned int _n_OP_vars;
  unsigned int _OP_number;

  std::vector<VariableValue *> _coupled_OP_vars;

  // Real _Galphamix;
//  Real _Gdeltamix;
};

#endif //ACCOUPLEDCALPHAD_H
