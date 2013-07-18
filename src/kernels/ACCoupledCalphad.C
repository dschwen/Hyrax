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
#include <cmath>

template<>
InputParameters validParams<ACCoupledCalphad>()
{
  InputParameters params = validParams<ACBulk>();
  params.addRequiredCoupledVar("coupled_CH_var", "the concentration to be coupled to the AC equation");
  params.addRequiredParam<int>("n_OP_vars", "# of coupled OP variables");
  params.addRequiredCoupledVar("OP_var_names", "Array of coupled OP variable names");
  params.addRequiredParam<int>("OP_number","# of the order parameter for this kernel, starting from 1");
  params.addRequiredParam<Real>("temperature", "Simulation temperature");
  params.addParam<Real>("gas_constant", 8.3144621, "Universal gas constant");
  //universal gas constant supplied here in J/mol-K
  params.addRequiredParam<Real>("well_height", "Free energy well height");

  return params;
}

ACCoupledCalphad::ACCoupledCalphad(const std::string & name, InputParameters parameters)
    : ACBulk(name, parameters),
      _G_hcp_Zr(getMaterialProperty<Real>("G_hcp_Zr")),
      _G_hcp_ZrH(getMaterialProperty<Real>("G_hcp_ZrH")),
      _G_fcc_Zr(getMaterialProperty<Real>("G_fcc_Zr")),
      _G_fcc_ZrH2(getMaterialProperty<Real>("G_fcc_ZrH2")),
      _G_H2(getMaterialProperty<Real>("G_H2")),
      _L0(getMaterialProperty<Real>("L0")),
      _L1(getMaterialProperty<Real>("L1")),

      _molarVol_alpha_Zr(getMaterialProperty<Real>("molar_volume_alpha_Zr")),
      _molarVol_delta_ZrH2(getMaterialProperty<Real>("molar_volume_delta_ZrH2")),

      _T(getParam<Real>("temperature")),
      _R(getParam<Real>("gas_constant")),
      _w(getParam<Real>("well_height")),

      _X(coupledValue("coupled_CH_var")),
      _n_OP_vars(getParam<int>("n_OP_vars")),
      _OP_number(getParam<int>("OP_number"))//,

      //_Galphamix(getMaterialProperty<Real>("Galphamix")),
      //_Gdeltamix(getMaterialProperty<Real>("Gdeltamix"))
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
  Real square_sum, quad_sum, square_mult;
  square_sum = quad_sum = 0.0;

  if (_n_OP_vars == 1)
    square_mult = 0.0;
  else
    square_mult = 1.0;

  Real dgdn, dHeavisidedn, d2gdn2, d2Heavisidedn2;

  // compute the coupled OP terms
  for(unsigned int i=0; i<_n_OP_vars; i++)
  {
    if(i != _OP_number-1)
    {
      Real OP;
      OP = (*_coupled_OP_vars[i])[_qp];
      //square_sum += ((*_coupled_OP_vars[i])[_qp])*((*_coupled_OP_vars[i])[_qp]);

      //quad_sum += ((*_coupled_OP_vars[i])[_qp])*((*_coupled_OP_vars[i])[_qp])*((*_coupled_OP_vars[i])[_qp])*((*_coupled_OP_vars[i])[_qp]);

      //square_mult *= ((*_coupled_OP_vars[i])[_qp])*((*_coupled_OP_vars[i])[_qp]);
      if (OP < 0.0)
      {
        square_sum += 0.0;
        quad_sum += 0.0;
        square_mult *= 0.0;
      }
      else if (OP > 1.0)
      {
        square_sum += 1.0;
        quad_sum += 1.0;
        square_mult *= 1.0;
      }
      else
      {
        square_sum += OP*OP;
        quad_sum += OP*OP*OP*OP;
        square_mult += OP*OP;
      }
    }
  }

  switch (type)
  {
  case Residual:

    dgdn = computeDBarrierDOP(square_sum, quad_sum, square_mult);

    dHeavisidedn = computeDHeavisideDOP();

    return (-computeGalphamix() + computeGdeltamix() )*dHeavisidedn + _w*dgdn;

  case Jacobian:

    d2gdn2 = computeD2BarrierDOP2(square_sum, quad_sum, square_mult);

   d2Heavisidedn2 = computeD2HeavisideDOP2();

    // return 0;
    // return _phi[_j][_qp]*(6-12*_u[_qp])*(-1*computeGalphamix() + computeGdeltamix()) + _w*(2 - 12*_u[_qp] + 12*_u[_qp]*_u[_qp] + 2*(square_sum) + 2*(quad_sum) + 12*_u[_qp]*_u[_qp]*(square_sum) + 2*square_mult)*_phi[_j][_qp];

   return  (-computeGalphamix() + computeGdeltamix())*d2Heavisidedn2*_phi[_j][_qp] + _w*d2gdn2*_phi[_j][_qp];
   // return (-1*_Galphamix[_qp] + _Gdeltamix[_qp])*d2Heavisidedn2*_phi[_j][_qp] + _w*d2gdn2*_phi[_j][_qp];
  }
  mooseError("Invalid type passed in");
}


Real
ACCoupledCalphad::computeDHeavisideDOP()
{
  if( _u[_qp] < 0.0 || _u[_qp] > 1.0 )
    return 0.0;
  else
    return  6.*_u[_qp]*(1. - _u[_qp]);
}


Real
ACCoupledCalphad::computeD2HeavisideDOP2()
{
  if( _u[_qp] < 0.0 || _u[_qp] > 1.0 )
    return 0.0;
  else
    return 6. - 12.*_u[_qp];
}


Real
ACCoupledCalphad::computeDBarrierDOP(Real & SS, Real & QS, Real & SM)
{
  if( _u[_qp] < 0.0 || _u[_qp] > 1.0 )
    return 0.;
  else
    return 2.*_u[_qp] - 6.*_u[_qp]*_u[_qp] + 4.*_u[_qp]*_u[_qp]*_u[_qp]
      + 2.*_u[_qp]*SS + 2.*_u[_qp]*QS + 4.*_u[_qp]*_u[_qp]*_u[_qp]*SS + 2.*_u[_qp]*SM;
}


Real
ACCoupledCalphad::computeD2BarrierDOP2(Real & SS, Real & QS, Real & SM)
{
  if( _u[_qp] < 0.0 || _u[_qp] > 1.0 )
    return 0.;
  else
    return 2.*(1. + QS + SM) + 12.*_u[_qp]*(-1. + _u[_qp] + _u[_qp]*SS);
}

Real
ACCoupledCalphad::computeGalphamix()
{
  Real reference, ideal;

  //do some checking so that equations are only calculated over the valid region
  //here, 0<atomic fraction H <= 0.5

  if(_X[_qp] > 0.5)
    return 0.0;
  else if(_X[_qp] < 0.)
    return 1e20;

  reference = (1. - 2.*_X[_qp])*_G_hcp_Zr[_qp] + _X[_qp]*_G_hcp_ZrH[_qp];

  ideal = _R*_T*( (1. - 2.*_X[_qp])*std::log( (2.*_X[_qp] - 1.)/(_X[_qp] - 1.))
                  +_X[_qp]*std::log( _X[_qp]/(1. - _X[_qp])) );

  return  (reference + ideal)/_molarVol_alpha_Zr[_qp];
}

Real
ACCoupledCalphad::computeGdeltamix()
{
  Real reference, ideal, excess;

  //do some checking so that equations are only calculated over the valid region
  //here, 0<atomic fraction H <= 2/3

  if(_X[_qp] > 2./3.)
    return 1e20;
  else if(_X[_qp] < 0.)
    return 0.0;

  reference = (0.5)*( (2. - 3.*_X[_qp])*_G_fcc_Zr[_qp] + _X[_qp]*_G_fcc_ZrH2[_qp] );

  ideal = _R*_T*( (2. - 3.*_X[_qp])*std::log( (2. - 3.*_X[_qp])/(2. - 2.*_X[_qp]))
                  +_X[_qp]*std::log( _X[_qp]/(2. - 2.*_X[_qp])) );

  excess = ((3.*_X[_qp] - 2.)*_X[_qp])/(4.*std::pow(_X[_qp] - 1., 2.));
  excess *= ( (_X[_qp] - 1.)*_L0[_qp] + (1. - 2.*_X[_qp])*_L1[_qp] );

  return (reference + ideal + excess)/_molarVol_delta_ZrH2[_qp];
}
