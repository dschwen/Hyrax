/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  13 February 2014
*
*************************************************************************/

#include "PrecipitateMatrixMisfitMaterial.h"

template<>
InputParameters validParams<PrecipitateMatrixMisfitMaterial>()
{
  InputParameters params = validParams<LinearSingleCrystalPrecipitateMaterial>();
  // matrix material into C_ijkl
  // precipitate material into C_precipitate
  // precipitate misfit into e_precipitate
  // number of precip variants into n_variants
  // order parameter coupled variables into variable_names
  // scaling factor for energy nondimensionalization into scaling_factor
  // temperature dependence of precip misfit coeffs into misfit_temperature_coeffs

  params.addRequiredParam<std::vector<Real> >("e_matrix","Eigenstrain tensor for solute in matrix: e11, e22, e33, e23, e13, e12");
  params.addRequiredCoupledVar("solute_name","coupled variable name of solute in matrix");

  return params;
}

PrecipitateMatrixMisfitMaterial::PrecipitateMatrixMisfitMaterial(const std::string & name,
                                                                 InputParameters parameters) :
    LinearSingleCrystalPrecipitateMaterial(name, parameters),
    _eigenstrain_matrix_vector(getParam<std::vector<Real> >("e_matrix")),
    _eigenstrain_matrix(),

    _eigenstrain_matrix_MP(declareProperty<RankTwoTensor>("eigenstrain_matrix_MP")),
    _dn_eigenstrain_matrix_MP(declareProperty<std::vector<RankTwoTensor> >("dn_eigenstrain_matrix_MP")),
    _dc_eigenstrain_matrix_MP(declareProperty<RankTwoTensor>("dc_eigenstrain_matrix_MP")),

    _dn_elasticity_tensor(declareProperty<std::vector<ElasticityTensorR4> >("dn_elasticity_tensor")),
    _dc_elasticity_tensor(declareProperty<ElasticityTensorR4>("dc_elasticity_tensor")),
    _dndn_elasticity_tensor(declareProperty<std::vector<ElasticityTensorR4> >("dndn_elasticity_tensor")),
    _dcdc_elasticity_tensor(declareProperty<ElasticityTensorR4>("dcdc_elasticity_tensor")),
    _dcdn_elasticity_tensor(declareProperty<std::vector<ElasticityTensorR4> >("dcdn_elasticity_tensor")),

    _dn_misfit_strain(declareProperty<std::vector<RankTwoTensor> >("dn_misfit_strain")),
    _dc_misfit_strain(declareProperty<RankTwoTensor>("dc_misfit_strain")),
    _dndn_misfit_strain(declareProperty<std::vector<RankTwoTensor> >("dndn_misfit_strain")),
    _dcdc_misfit_strain(declareProperty<RankTwoTensor>("dcdc_misfit_strain")),
    _dcdn_misfit_strain(declareProperty<std::vector<RankTwoTensor> >("dcdn_misfit_strain")),

    _solute(coupledValue("solute_name"))
{
  // fill in the original tensors.  Don't touch after this!
  _eigenstrain_matrix.fillFromInputVector(_eigenstrain_matrix_vector);
}

void
PrecipitateMatrixMisfitMaterial::computeProperties()
{
  for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
  {
    _eigenstrains_MP[_qp].resize(_n_variants);
    _d_eigenstrains_MP[_qp].resize(_n_variants);
    _dn_eigenstrain_matrix_MP[_qp].resize(_n_variants);

    _precipitate_eigenstrain[_qp].resize(_n_variants);

    _dn_elasticity_tensor[_qp].resize(_n_variants);
    _dndn_elasticity_tensor[_qp].resize(_n_variants);
    _dcdn_elasticity_tensor[_qp].resize(_n_variants);

    _dn_misfit_strain[_qp].resize(_n_variants);
    _dcdn_misfit_strain[_qp].resize(_n_variants);
    _dndn_misfit_strain[_qp].resize(_n_variants);
    _dcdn_misfit_strain[_qp].resize(_n_variants);

    computeQpElasticityTensor();
    computeQpEigenstrain();
    computeQpElasticStrain();
    computeQpElasticStress();
  }
}

void
PrecipitateMatrixMisfitMaterial::computeQpElasticityTensor()
{
  //going to need to put temperature dependence in here

  Real inverse = 1/_scaling_factor;
  ElasticityTensorR4 zeros;
  zeros.zero();

  _Cijkl_MP[_qp] = _Cijkl/inverse;
  _Cijkl_precipitates_MP[_qp] = _Cijkl_precipitate/inverse;

  Real sum_OP = 0;
  for (unsigned int i=0; i<_n_variants; i++)
  {
    sum_OP += (*_OP[i])[_qp]*(*_OP[i])[_qp];
    //(_dn_elasticity_tensor[_qp])[i] = ( (_Cijkl_precipitate - _Cijkl)*2*(*_OP[i])[_qp] )/inverse;
    (_dn_elasticity_tensor[_qp])[i] = (_Cijkl_precipitate - _Cijkl)*computeDHeaviside(i)/inverse;
    //(_dndn_elasticity_tensor[_qp])[i] = ( (_Cijkl_precipitate - _Cijkl)*2 )/inverse;
    (_dndn_elasticity_tensor[_qp])[i] = ( (_Cijkl_precipitate - _Cijkl)*computeD2Heaviside(i) )/inverse;
    (_dcdn_elasticity_tensor[_qp])[i] = zeros;
  }

  //_elasticity_tensor[_qp] = (_Cijkl + (_Cijkl_precipitate - _Cijkl)*sum_OP)/inverse;
  _elasticity_tensor[_qp] = (_Cijkl + (_Cijkl_precipitate - _Cijkl)*computeHeaviside())/inverse;
  _dc_elasticity_tensor[_qp] = zeros;
  _dcdc_elasticity_tensor[_qp] = zeros;

  //not sure if I actually got the Jacobian multiplier right here.
  //_Jacobian_mult[_qp] = (_Cijkl + (_Cijkl_precipitate - _Cijkl)*sum_OP)/inverse;
  _Jacobian_mult[_qp] = (_Cijkl + (_Cijkl_precipitate - _Cijkl)*computeHeaviside())/inverse;
}

void
PrecipitateMatrixMisfitMaterial::computeQpEigenstrain()
{
  computeQpPrecipitateEigenstrain();
  computeQpMatrixEigenstrain();
}

void
PrecipitateMatrixMisfitMaterial::computeQpPrecipitateEigenstrain()
{
  LinearSingleCrystalPrecipitateMaterial::computeQpEigenstrain();
  //might want to tweak the derivative values...
}

void
PrecipitateMatrixMisfitMaterial::computeQpMatrixEigenstrain()
{
  Real interpolation_value(0);

  for (unsigned int i=0; i<_n_variants; i++)
  {
    interpolation_value =+ (*_OP[i])[_qp]*(*_OP[i])[_qp];
    (_dn_eigenstrain_matrix_MP[_qp])[i] = _eigenstrain_matrix*2*_solute[_qp]*(*_OP[i])[_qp];
  }

  _eigenstrain_matrix_MP[_qp] = _eigenstrain_matrix*_solute[_qp]*(1 - interpolation_value);
  _dc_eigenstrain_matrix_MP[_qp] = _eigenstrain_matrix*(1 - interpolation_value);

  //miiight want some stuff here for derivative values...

}

void
PrecipitateMatrixMisfitMaterial::computeQpElasticStrain()
{
  // compute the elastic strain: e_el = e_local - e_misfit

  //compute the total strain
  computeQpStrain();

  //compute the overall misfit strain
  computeQpMisfitStrain();

  _elastic_strain[_qp] = _local_strain[_qp] - _misfit_strain[_qp];
}

void
PrecipitateMatrixMisfitMaterial::computeQpMisfitStrain()
{
  //sum up the misfit strains for the orientation variants
  RankTwoTensor sum_precipitate_strains;
  RankTwoTensor zeros;
  sum_precipitate_strains.zero();
  zeros.zero();

  Real OP_sum = 0;

  for(unsigned int i=0; i<_n_variants; i++)
  {
    OP_sum += (*_OP[i])[_qp]*(*_OP[i])[_qp];
    sum_precipitate_strains += (_eigenstrains_MP[_qp])[i];
    (_dn_misfit_strain[_qp])[i] = (_d_eigenstrains_MP[_qp])[i] - (_dn_eigenstrain_matrix_MP[_qp])[i];

    (_dcdn_misfit_strain[_qp])[i] = ( (_precipitate_eigenstrain[_qp])[i] -_eigenstrain_matrix*_solute[_qp] )*-2;
    (_dndn_misfit_strain[_qp])[i] = _eigenstrain_matrix*2*(*_OP[i])[_qp];
  }

  _misfit_strain[_qp] = sum_precipitate_strains + _eigenstrain_matrix_MP[_qp];
  _dc_misfit_strain[_qp] = _eigenstrain_matrix*(1 - OP_sum);
  _dcdc_misfit_strain[_qp] = zeros;
}


//void
//
//  computeQpStress();
//}

//void
//PrecipitateMatrixMisfitMaterial::computeQpStrain()
//{
//  //strain = (grad_disp + grad_disp^T)/2
//  RankTwoTensor grad_tensor(_grad_disp_x[_qp],_grad_disp_y[_qp],_grad_disp_z[_qp]);

//  _local_strain[_qp] = (grad_tensor + grad_tensor.transpose())/2.0;
//}

//void
//PrecipitateMatrixMisfitMaterial::computeQpStress()
//{
//  // stress = C * e
//  _stress[_qp] = _elasticity_tensor[_qp]*_elastic_strain[_qp];
//}


Real
PrecipitateMatrixMisfitMaterial::computeHeaviside()
{
  Real heaviside_first(0);
  Real heaviside_second(0);

  //may need to put some checking in here so that OP fixed between 0 and 1
  for(unsigned int i=0; i<_n_variants; i++)
  {
    heaviside_first += std::pow((*_OP[i])[_qp], 2);
    heaviside_second += std::pow((*_OP[i])[_qp], 3);
  }

  return 3*heaviside_first - 2*heaviside_second;
}

Real
PrecipitateMatrixMisfitMaterial::computeDHeaviside(unsigned int i)
{
  return 6*(*_OP[i])[_qp]*(1 - (*_OP[i])[_qp]);
}

Real
PrecipitateMatrixMisfitMaterial::computeD2Heaviside(unsigned int i)
{
  return 6*(1 - (*_OP[i])[_qp]);
}
