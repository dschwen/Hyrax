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

  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

private:
  MaterialProperty<Real> & _W;                            //Well height
  MaterialProperty<Real> & _Omega;                        //Molar volume
  MaterialProperty<Real> & _G_alpha;                      //Gmix_alpha
  MaterialProperty<Real> & _G_delta;                      //Gmix_delta
  MaterialProperty<Real> & _dGalpha_dc;
  MaterialProperty<Real> & _dGdelta_dc;

  unsigned int _n_OP_vars;
  unsigned int _OP_number;

  std::vector<unsigned int> _n_var;
  std::vector<VariableValue *> _coupled_OP_vars;

  Real _scaling_factor;

  unsigned int _c_var;
  unsigned int _w_var;
  unsigned int _T_var;

  VariableValue & _c;
  VariableValue & _w;
  VariableValue & _T;

};

#endif //ACCOUPLEDCALPHAD_H
