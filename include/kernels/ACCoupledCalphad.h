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

private:
  MaterialProperty<Real> & _W;                            //Well height
  MaterialProperty<Real> & _Omega;                        //Molar volume
  MaterialProperty<Real> & _G_alpha;                      //Gmix_alpha
  MaterialProperty<Real> & _G_delta;                      //Gmix_delta

  unsigned int _n_OP_vars;
  unsigned int _OP_number;

  std::vector<VariableValue *> _coupled_OP_vars;

};

#endif //ACCOUPLEDCALPHAD_H
