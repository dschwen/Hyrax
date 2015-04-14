/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*
*  19 November 2013
*
*************************************************************************/

#include "CHCoupledCalphadSplit.h"

template<>
InputParameters validParams<CHCoupledCalphadSplit>()
{
  InputParameters params = validParams<SplitCHCRes>();

  //params.addRequiredCoupledVar("coupled_OP_var", "The order parameter coupled to the CH eqn");
  params.addParam<Real>("scaling_factor", 1, "free energy scaling factor for nondimensionalization");

  // params.addRequiredCoupledVar("n", "structural order parameter");
  params.addRequiredCoupledVar("T", "temperature");

  params.addRequiredParam<int>("n_OP_vars", "# of coupled OP variables");
  params.addRequiredCoupledVar("OP_var_names", "Array of coupled OP variable names");

  return params;
}

CHCoupledCalphadSplit::CHCoupledCalphadSplit(const std::string & name, InputParameters parameters) :
    SplitCHCRes(name, parameters),
    _W(getMaterialProperty<Real>("well_height")),
    _Omega(getMaterialProperty<Real>("molar_volume")),
    _dGalpha_dc(getMaterialProperty<Real>("dGAB1CD1_dc")),
    _d2Galpha_dc2(getMaterialProperty<Real>("d2GAB1CD1_dc2")),
    _dGdelta_dc(getMaterialProperty<Real>("dGAB1CD2_dc")),
    _d2Gdelta_dc2(getMaterialProperty<Real>("d2GAB1CD2_dc2")),
    // _OP(coupledValue("coupled_OP_var")),
    //_OP_var(coupled("coupled_OP_var")),
    _scaling_factor(getParam<Real>("scaling_factor")),
    // _n_var(coupled("n")),
    _w_var(coupled("w")),
    _T_var(coupled("T")),
    //_n(coupledValue("n")),
    _w(coupledValue("w")),
    _T(coupledValue("T")),
    _n_OP_vars(getParam<int>("n_OP_vars"))
{
  if(_n_OP_vars != coupledComponents("OP_var_names"))
    mooseError("Please match the number of orientation variants to coupled OPs (ACCoupledCalphad).");

  _n_var.resize(_n_OP_vars);
  _OP.resize(_n_OP_vars);

  for (unsigned int i=0; i<_n_OP_vars; i++)
  {
    _n_var[i] = coupled("OP_var_names", i);
    _OP[i] = &coupledValue("OP_var_names", i);
  }
}

Real
CHCoupledCalphadSplit::computeDFDC(PFFunctionType type)
{
  Real Heaviside, dHeaviside;

  Heaviside = computeHeaviside();

  switch (type)
  {
  case Residual:
    return _scaling_factor*
      ( ( (1 - Heaviside)*_dGalpha_dc[_qp] + Heaviside*_dGdelta_dc[_qp] )/_Omega[_qp]);

  case Jacobian:
    return _scaling_factor*
      ( _phi[_j][_qp]*((1 - Heaviside)*_d2Galpha_dc2[_qp] + Heaviside*_d2Gdelta_dc2[_qp] )/_Omega[_qp]);

  }
  mooseError("invalid type passed in");
}

Real
CHCoupledCalphadSplit::computeQpOffDiagJacobian(unsigned int jvar)
{
  for (unsigned int i=0; i<_n_OP_vars; i++)
  {
    if (jvar == _n_var[i])
        return _scaling_factor*_phi[_j][_qp]*_test[_i][_qp]*
          ( computeDHeaviside(i)*(_dGdelta_dc[_qp] - _dGalpha_dc[_qp]) )/_Omega[_qp];
  }



//  if (jvar == _n_var)
  //  return _scaling_factor*_phi[_j][_qp]*_test[_i][_qp]*
  //  ( computeDHeaviside()*(_dGdelta_dc[_qp] - _dGalpha_dc[_qp]) )/_Omega[_qp];

 if (jvar == _w_var)
    return -_phi[_j][_qp] * _test[_i][_qp];

  else if (jvar == _T_var)
    return 0;

  else
    mooseError("Screwed up CHCoupledCalphadSplit::computeQpOffDiagJacobian()");
}

Real
CHCoupledCalphadSplit::computeHeaviside()
{
  Real heaviside_first(0);
  Real heaviside_second(0);

  //may need to put some checking in here so that OP fixed between 0 and 1
  for(unsigned int i=0; i<_n_OP_vars; i++)
  {
    heaviside_first += std::pow((*_OP[i])[_qp], 2);
    heaviside_second += std::pow((*_OP[i])[_qp], 3);
  }

  // return 3*heaviside_first - 2*heaviside_second;
  //testing to see if n works, instead of heaviside
  return (*_OP[i])[_qp];
}

Real
CHCoupledCalphadSplit::computeDHeaviside(unsigned int i)
{
  //return 6*(*_OP[i])[_qp]*(1 - (*_OP[i])[_qp]);
  //testing to see if n works, instead of heaviside
  return 0;
}
