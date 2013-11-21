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
  // params.addRequiredParam<Real>("Well_height", "well height of the double-well functional");
  params.addRequiredCoupledVar("coupled_OP_var", "The order parameter coupled to the CH eqn");

  return params;
}

CHCoupledCalphadSplit::CHCoupledCalphadSplit(const std::string & name, InputParameters parameters):
    SplitCHCRes(name, parameters),
    _W(getMaterialProperty<Real>("well_height")),
    _Omega(getMaterialProperty<Real>("molar_volume")),
    _dGalpha_dc(getMaterialProperty<Real>("dGAB1CD1_dc")),
    _d2Galpha_dc2(getMaterialProperty<Real>("d2GAB1CD1_dc2")),
    _dGdelta_dc(getMaterialProperty<Real>("dGAB1CD2_dc")),
    _d2Gdelta_dc2(getMaterialProperty<Real>("d2GAB1CD2_dc2")),
    _OP(coupledValue("coupled_OP_var"))
{
}


Real
CHCoupledCalphadSplit::computeDFDC(PFFunctionType type)
{
  Real Heaviside, dHeaviside;

  Heaviside = computeHeaviside();
  // dHeaviside = computeDHeaviside();

  switch (type)
  {
  case Residual:
    return ( (1 - Heaviside)*_dGalpha_dc[_qp] + Heaviside*_dGdelta_dc[_qp] )/_Omega[_qp];


  case Jacobian:

    return _phi[_j][_qp]*((1 - Heaviside)*_d2Galpha_dc2[_qp] + Heaviside*_d2Gdelta_dc2[_qp] )/_Omega[_qp];

  }
  mooseError("invalid type passed in");
}

Real
CHCoupledCalphadSplit::computeHeaviside()
{
  Real heaviside_first(0);
  Real heaviside_second(0);

  //may need to put some checking in here so that OP fixed between 0 and 1
  //for(unsigned int i=0; i<_n_OP_variables; i++)
  //{
    heaviside_first += std::pow(_OP[_qp], 2);
    heaviside_second += std::pow(_OP[_qp], 3);
    //}

  return 3*heaviside_first - 2*heaviside_second;
}
