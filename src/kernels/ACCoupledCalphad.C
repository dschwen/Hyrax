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
  params.addParam<Real>("scaling_factor", 1, "free energy scaling factor for nondimensionalization");

  params.addRequiredCoupledVar("c", "concentration");
  params.addRequiredCoupledVar("w", "chemical potential");
  params.addRequiredCoupledVar("T", "temperature");

  return params;
}

ACCoupledCalphad::ACCoupledCalphad(const std::string & name, InputParameters parameters) :
    ACBulk(name, parameters),
    _W(getMaterialProperty<Real>("well_height")),
    _Omega(getMaterialProperty<Real>("molar_volume")),
    _G_alpha(getMaterialProperty<Real>("G_AB1CD1")),
    _G_delta(getMaterialProperty<Real>("G_AB1CD2")),
    _dGalpha_dc(getMaterialProperty<Real>("dGAB1CD1_dc")),
    _dGdelta_dc(getMaterialProperty<Real>("dGAB1CD2_dc")),

    _n_OP_vars(getParam<int>("n_OP_vars")),
    _OP_number(getParam<int>("OP_number")),
    _scaling_factor(getParam<Real>("scaling_factor")),

    _c_var(coupled("c")),
    _w_var(coupled("w")),
    _T_var(coupled("T")),

    _c(coupledValue("c")),
    _w(coupledValue("w")),
    _T(coupledValue("T"))
{
  // Create a vector of the coupled OP variables and set = 0 the one that the kernel
  // is operating on
  if(_n_OP_vars != coupledComponents("OP_var_names"))
    mooseError("Please match the number of orientation variants to coupled OPs (ACCoupledCalphad).");

  _n_var.resize(_n_OP_vars);
  _coupled_OP_vars.resize(_n_OP_vars);

  for(unsigned int i=0; i< _n_OP_vars; i++)
  {
    if(i == _OP_number-1)
      _coupled_OP_vars[i] = NULL;
    else
      _coupled_OP_vars[i] = &coupledValue("OP_var_names", i);

    if (i != _OP_number-1)
      _n_var[i] = coupled("OP_var_names", i);
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

  //compute the coupled OP terms
  for(unsigned int i=0; i<_n_OP_vars; i++)
  {
    if(i != _OP_number-1)
    {
      Real OP;
      OP = (*_coupled_OP_vars[i])[_qp];

      square_sum += OP*OP;
      quad_sum += OP*OP*OP*OP;
      square_mult *= OP*OP;
    }
  }

  // _console<<"OP number = "<<_OP_number<<std::endl;
  // _console<<"square_sum = "<<square_sum<<std::endl;
  // _console<<"square_mult = "<<square_mult<<std::endl;
  // _console<<"quad_sum = "<<quad_sum<<std::endl<<std::endl;

  Real dgdn, dHeavisidedn, d2gdn2, d2Heavisidedn2;
  switch (type)
  {
  case Residual:

    dgdn = computeDBarrierDOP(square_sum, quad_sum, square_mult);

    dHeavisidedn = computeDHeavisideDOP();

    return _scaling_factor*( ( (_G_delta[_qp] - _G_alpha[_qp])*dHeavisidedn + _W[_qp]*dgdn ) / _Omega[_qp] );

  case Jacobian:

    d2gdn2 = computeD2BarrierDOP2(square_sum, quad_sum, square_mult);

    d2Heavisidedn2 = computeD2HeavisideDOP2();

    return  _scaling_factor*( ( ( _G_delta[_qp] - _G_alpha[_qp])*d2Heavisidedn2*_phi[_j][_qp]
                 + _W[_qp]*d2gdn2*_phi[_j][_qp] ) / _Omega[_qp] );
  }
  mooseError("Invalid type passed in");
}

Real
ACCoupledCalphad::computeQpOffDiagJacobian(unsigned int jvar)
{
  /* for (unsigned int i=0; i< _n_OP_vars; i++)
    {
    if (i != _OP_number)
    {
      if (jvar == _n_var[i])
      {
        //Real val = (*_coupled_OP_vars[i])[_qp];
        //Real s = (*_coupled_OP_vars[_n_OP_vars - (i+_OP_number)%_n_OP_vars -1])[_qp];

        //Real ans = 2*val*val*(2*_u[_qp] + 2*val*val*_u[_qp] + 4*_u[_qp]*_u[_qp]*_u[_qp] + 2*_u[_qp]*s*s );
        //ans /= _Omega[_qp];
        //ans *= _test[_i][_qp]*_phi[_j][_qp];
        //return ans;
        return 0;
      }

    }
    }*/

  if (jvar == _c_var)
    return _L[_qp]*_test[_i][_qp]*_phi[_j][_qp]*
      (_scaling_factor*( computeDHeavisideDOP()*(_dGdelta_dc[_qp] - _dGalpha_dc[_qp]) )/_Omega[_qp] );

  else if (jvar == _w_var)
    return 0;

  else if (jvar == _T_var)
    //return _L[_qp]*_test[_i][_qp]*_phi[_j][_qp]*(d/dT of df/dn);
    return 0;

  else
    return 0;
    // mooseError("Screwed up ACCoupledCalphad::computeQpOffDiagJacobian()");
}


Real
ACCoupledCalphad::computeDHeavisideDOP()
{
  //return  6*_u[_qp]*(1 - _u[_qp]);

  Real OP =_u[_qp];
  if(OP < 0) OP = 0;
  if(OP > 1) OP = 1;

  return 6*OP*(1-OP);

  //testing h as 7th order polynomial
  // return -56.59203*7*OP*OP*OP*OP*OP*OP + 198.28747*6*OP*OP*OP*OP*OP -277.2044*5*OP*OP*OP*OP + 196.95257*4*OP*OP*OP -74.5349*3*OP*OP + 14.1276*2*OP - 0.036161;

  //testing h as 5th order polynomial

  //return -10.957*5*OP*OP*OP*OP + 28.1*4*OP*OP*OP - 25.5071*3*OP*OP - 9.616*2*OP - 0.25672;

}


Real
ACCoupledCalphad::computeD2HeavisideDOP2()
{
  //return 6 - 12*_u[_qp];

  Real OP = _u[_qp];
  if(OP < 0) return 0;
  if(OP > 1) return 0;

  return 6 - 12*OP;


//testing h as 7th order polynomial
   //return -56.59203*7*6*OP*OP*OP*OP*OP + 198.28747*6*5*OP*OP*OP*OP -277.2044*5*4*OP*OP*OP  + 196.95257*4*3*OP*OP -74.5349*3*2*OP + 14.1276*2;

  //testing h as 5th order polynomial

//  return -10.957*5*4*OP*OP*OP + 28.1*4*3*OP*OP - 25.5071*3*2*OP - 9.616*2;

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
