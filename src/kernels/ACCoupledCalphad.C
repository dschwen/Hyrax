/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  13 June 2013
*
*************************************************************************/

#include "ACCoupledCalphad.h"

template<>
InputParameters validParams<ACCoupledCalphad>()
{
  InputParameters params = validParams<ACBulk>();
  params.addRequiredParam<int>("n_OP_vars", "# of coupled OP variables");
  params.addRequiredCoupledVar("OP_var_names", "Array of coupled OP variable names");
  params.addRequiredParam<int>("OP_number","# of the order parameter for this kernel, starting from 1");

  return params;
}

ACCoupledCalphad::ACCoupledCalphad(const std::string & name, InputParameters parameters)
    : ACBulk(name, parameters),
      _W(getMaterialProperty<Real>("well_height")),
      _Omega(getMaterialProperty<Real>("molar_volume")),
      _G_alpha(getMaterialProperty<Real>("G_AB1CD1")),
      _G_delta(getMaterialProperty<Real>("G_AB1CD2")),

      _n_OP_vars(getParam<int>("n_OP_vars")),
      _OP_number(getParam<int>("OP_number"))
{
  // Create a vector of the coupled OP variables and set = 0 the one that the kernel
  // is operating on
  if(_n_OP_vars != coupledComponents("OP_var_names"))
    mooseError("Please match the number of orientation variants to coupled OPs (ACCoupledCalphad).");

  _coupled_OP_vars.resize(_n_OP_vars);

  for(unsigned int i=0; i< _n_OP_vars; i++)
  {
    if(i == _OP_number-1)
      _coupled_OP_vars[i] = NULL;
    else
      _coupled_OP_vars[i] = &coupledValue("OP_var_names", i);
  }
}

Real
ACCoupledCalphad::computeDFDOP(PFFunctionType type)
{
  //std::cout<<"in ACCoupledCalphad computeDFDOP"<<std::endl;
  Real square_sum, quad_sum, square_mult;
  square_sum = quad_sum = 0.0;

  if (_n_OP_vars == 1)
    square_mult = 0.0;
  else
    square_mult = 1.0;

  //compute the coupled OP terms
  for(unsigned int i=0; i<_n_OP_vars; i++)
  {
    if(i != _OP_number-1)
    {
      Real OP;
      OP = (*_coupled_OP_vars[i])[_qp];

      square_sum += OP*OP;
      quad_sum += OP*OP*OP*OP;
      square_mult += OP*OP;
    }
  }

  Real dgdn, dHeavisidedn, d2gdn2, d2Heavisidedn2;
  switch (type)
  {
  case Residual:

    dgdn = computeDBarrierDOP(square_sum, quad_sum, square_mult);

    dHeavisidedn = computeDHeavisideDOP();

    return ( (_G_delta[_qp] - _G_alpha[_qp])*dHeavisidedn + _W[_qp]*dgdn ) / _Omega[_qp];

  case Jacobian:

    d2gdn2 = computeD2BarrierDOP2(square_sum, quad_sum, square_mult);

    d2Heavisidedn2 = computeD2HeavisideDOP2();

    return  ( ( _G_delta[_qp] - _G_alpha[_qp])*d2Heavisidedn2*_phi[_j][_qp]
              + _W[_qp]*d2gdn2*_phi[_j][_qp] ) / _Omega[_qp];
  }
  mooseError("Invalid type passed in");
}

Real
ACCoupledCalphad::computeDHeavisideDOP()
{
  return  6*_u[_qp]*(1 - _u[_qp]);
}


Real
ACCoupledCalphad::computeD2HeavisideDOP2()
{
  return 6 - 12*_u[_qp];
}

Real
ACCoupledCalphad::computeDBarrierDOP(Real & SS, Real & QS, Real & SM)
{
  return 2*_u[_qp] - 6*_u[_qp]*_u[_qp] + 4*_u[_qp]*_u[_qp]*_u[_qp]
      + 2*_u[_qp]*SS + 2*_u[_qp]*QS + 4*_u[_qp]*_u[_qp]*_u[_qp]*SS + 2*_u[_qp]*SM;
}

Real
ACCoupledCalphad::computeD2BarrierDOP2(Real & SS, Real & QS, Real & SM)
{
  return 2 - 12*_u[_qp] + 12*_u[_qp]*_u[_qp] + 2*SS + 2*QS + 12*_u[_qp]*_u[_qp]*SS + 2*SM;
}
