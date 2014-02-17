/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  14 February 2014
*
*************************************************************************/

#include "ACPrecipMatrixElasticity.h"

template<>
InputParameters validParams<ACPrecipMatrixElasticity>()
{
  InputParameters params = validParams<ACBulk>();
  params.addRequiredParam<int>("OP_number","# of the order parameter for this kernel, starting from 1");
 params.addParam<Real>("scaling_factor", 1, "elastic energy scaling factor for nondimensionalization");

  return params;
}

ACPrecipMatrixElasticity::ACPrecipMatrixElasticity(const std::string & name, InputParameters parameters)
    : ACBulk(name, parameters),
      _elasticity_tensor(getMaterialProperty<ElasticityTensorR4>("elasticity_tensor")),
      _dn_elasticity_tensor(getMaterialProperty<std::vector<ElasticityTensorR4> >("dn_elasticity_tensor")),
      _elastic_strain(getMaterialProperty<RankTwoTensor>("elastic_strain")),
      _dn_misfit_strain(getMaterialProperty<std::vector<RankTwoTensor> >("dn_misfit_strain")),
      _OP_number(getParam<int>("OP_number")),
      _scaling_factor(getParam<Real>("scaling_factor"))
{
}

Real
ACPrecipMatrixElasticity::computeDFDOP(PFFunctionType type)
{
  RankTwoTensor a = _elasticity_tensor[_qp]*(_elastic_strain[_qp]);
  RankTwoTensor b = _elasticity_tensor[_qp]*(_dn_misfit_strain[_qp])[_OP_number-1]*(-1);
  RankTwoTensor c = (_dn_elasticity_tensor[_qp])[_OP_number-1]*_elastic_strain[_qp];

  Real first  = a.doubleContraction( (_dn_misfit_strain[_qp])[_OP_number-1]*(-1) );
  Real second = b.doubleContraction( _elastic_strain[_qp] );
  Real third  = c.doubleContraction( _elastic_strain[_qp] );

  switch(type)
  {
  case Residual:
    return _scaling_factor*0.5*(first + second + third);

  case Jacobian:
    return 0;
  }

  mooseError("invalid type passed in");
}

