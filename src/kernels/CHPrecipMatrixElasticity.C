/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  14 February 2014
*
*************************************************************************/

#include "CHPrecipMatrixElasticity.h"

template<>
InputParameters validParams<CHPrecipMatrixElasticity>()
{
  InputParameters params = validParams<SplitCHCRes>();

  params.addParam<Real>("scaling_factor", 1, "elastic energy scaling factor for nondimensionalization");

  return params;
}

CHPrecipMatrixElasticity::CHPrecipMatrixElasticity(const std::string & name, InputParameters parameters)
    : SplitCHCRes(name, parameters),
      _elasticity_tensor(getMaterialProperty<ElasticityTensorR4>("elasticity_tensor")),
      _elastic_strain(getMaterialProperty<RankTwoTensor>("elastic_strain")),
      _dc_misfit_strain(getMaterialProperty<RankTwoTensor>("dc_misfit_strain")),
      _scaling_factor(getParam<Real>("scaling_factor"))
{
}

Real
CHPrecipMatrixElasticity::computeDFDC(PFFunctionType type)
{
  RankTwoTensor a = _elasticity_tensor[_qp]*(_elastic_strain[_qp]);
  RankTwoTensor b = _elasticity_tensor[_qp]*_dc_misfit_strain[_qp]*(-1);

  Real first  = a.doubleContraction( (_dc_misfit_strain[_qp])*(-1) );
  Real second = b.doubleContraction( _elastic_strain[_qp] );

  switch(type)
  {
  case Residual:
    return _scaling_factor*0.5*(first + second);

  case Jacobian:
    return 0;
  }

  mooseError("invalid type passed in");
}


//really hope i don't need any off diagonal shit here
