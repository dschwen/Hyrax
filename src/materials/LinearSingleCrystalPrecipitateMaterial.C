/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  12 January 2012
*
*************************************************************************/

#include "LinearSingleCrystalPrecipitateMaterial.h"

/**
 * LinearSingleCrystalPrecipitateMaterial handles anisotropic,
 * single-crystal material elastic constants.  It handles a single
 * crystal of matrix with an arbitrary number of orientation variants of
 * a coherent precipitate.
 *
 * The matrix material goes into C_ijkl, while the precipitate info goes into
 * C_ijkl_precipitate
 */

template<>
InputParameters validParams<LinearSingleCrystalPrecipitateMaterial>()
{
  InputParameters params = validParams<TensorMechanicsMaterial>();
  // matrix material into C_ijkl
  params.addRequiredParam<std::vector<Real> >("C_precipitate", "Stiffness tensor for precipitate");
  params.addRequiredParam<std::vector<Real> >("e_precipitate","Eigenstrain tensor for precipitate: e11, e22, e33, e23, e13, e12");
  params.addRequiredParam<int>("n_variants","# of orientation variants for precipitate in single crystal");
  params.addRequiredCoupledVar("variable_names","Array of coupled variable names");

  return params;
}

LinearSingleCrystalPrecipitateMaterial::LinearSingleCrystalPrecipitateMaterial(const std::string & name, InputParameters parameters)
    : TensorMechanicsMaterial(name, parameters),
      _Cijkl_precipitate_vector(getParam<std::vector<Real> >("C_precipitate")),
      _eigenstrain_vector(getParam<std::vector<Real> >("e_precipitate")),

      _Cijkl_precipitate(),
      _eigenstrain(),
      _n_variants(getParam<int>("n_variants")),

      //  _Cijkl_precipitates_rotated(),
      _eigenstrains_rotated(),

      _local_strain(declareProperty<RankTwoTensor >("local_strain")),
      _misfit_strain(declareProperty<RankTwoTensor >("misfit_strain")),

      _eigenstrains_MP(declareProperty<std::vector<RankTwoTensor> >("eigenstrains_MP")),
      _Cijkl_MP(declareProperty<ElasticityTensorR4>("Cijkl_MP")),
      _Cijkl_precipitates_MP(declareProperty<ElasticityTensorR4>("Cijkl_precipitates_MP")),
      //_d_elasticity_tensor(declareProperty<std::vector<ElasticityTensorR4> >("d_elasticity_tensor")),
      _d_eigenstrains_MP(declareProperty<std::vector<RankTwoTensor> >("d_eigenstrains_MP")),
      _precipitate_eigenstrain(declareProperty<std::vector<RankTwoTensor> >("precipitate_eigenstrain"))
{
  // check to make sure the input file is all set up right
  if(_n_variants != coupledComponents("variable_names"))
    mooseError("Please match the number of orientation variants with coupled order parameters (LSXPM).");

  // size vectors appropriately
  _coupled_variables.resize(_n_variants);
  // _Cijkl_precipitates_rotated.resize(_n_variants);
  _eigenstrains_rotated.resize(_n_variants);

  // populate with data
  for(unsigned int i=0; i < _n_variants; i++)
    _coupled_variables[i] = &coupledValue("variable_names", i);

  _Cijkl_precipitate.fillFromInputVector(_Cijkl_precipitate_vector, _all_21);
  _eigenstrain.fillFromInputVector(_eigenstrain_vector);

  // fill in the first variant without rotation
  //_Cijkl_precipitates_rotated[0] = _Cijkl_precipitate;
  _eigenstrains_rotated[0] = _eigenstrain;

  // rotate all the things, in radians
  Real rotation_angle_base = 2.0*libMesh::pi/Real(_n_variants);
  Real rotation_angle = rotation_angle_base;

  //ElasticityTensorR4 Cijkl_to_rotate;
  //RankTwoTensor eigenstrain_to_rotate;

  //RotationTensor R(_Euler_angles);

  for(unsigned int i=1; i<_n_variants; i++)
  {
    // _Euler_angles(0) = rotation_angle;
    //R.update(_Euler_angles);

    // Cijkl_to_rotate = _Cijkl_precipitate;
    //Cijkl_to_rotate.rotate(R);
    //_Cijkl_precipitates_rotated[i] = Cijkl_to_rotate;
    //_Cijkl_precipitates_rotated[i] = _Cijkl_precipitate;

    //  eigenstrain_to_rotate = _eigenstrain;
    //eigenstrain_to_rotate.rotate(R);
    //_eigenstrains_rotated[i] = eigenstrain_to_rotate;
    _eigenstrains_rotated[i] = _eigenstrain.rotateXyPlane(rotation_angle);

    // increment the rotation angle for the next go-round
    rotation_angle = rotation_angle + rotation_angle_base;
  }
}

void
LinearSingleCrystalPrecipitateMaterial::computeProperties()
{
  for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
  {
    // resize all the material properties vectors.  Don't forget this.
    //  _Cijkl_precipitates_MP[_qp].resize(_n_variants);
    _eigenstrains_MP[_qp].resize(_n_variants);
    //_d_elasticity_tensor[_qp].resize(_n_variants);
    _d_eigenstrains_MP[_qp].resize(_n_variants);
    _precipitate_eigenstrain[_qp].resize(_n_variants);

    computeQpElasticityTensor();
    computeQpEigenstrain();
    computeQpElasticStrain();
    computeQpElasticStress();
  }
}

void
LinearSingleCrystalPrecipitateMaterial::computeQpElasticityTensor()
{
  /**
   * Not making homogeneous modulus approximation between matrix and precipitate.  Works only
   * if the OP goes between zero and one!  (Forcing it to max at one, right now).
   **/

  // Fill in the matrix stiffness material property
  _Cijkl_MP[_qp] = _Cijkl;

  // fill in the precipitate stiffnesses material property
  // for(unsigned int i(0); i < _n_variants; ++i)
  //  (_Cijkl_precipitates_MP[_qp])[i] = _Cijkl_precipitates_rotated[i];
  _Cijkl_precipitates_MP[_qp] = _Cijkl_precipitate;

  // Sum the order parameters and stiffnesses for the precipitates
  ElasticityTensorR4 sum_precipitate_tensors;

  sum_precipitate_tensors.zero();

  // interpolation for Cijkl dependence on order parameter
  Real interpolation_value;
  Real sum_interpolation_values(0.0);

  for(unsigned int i=0; i<_n_variants; i++)
  {
    // calculate the value of the interpolation function for an OP and this _qp
    interpolation_value = (*_coupled_variables[i])[_qp];
    // this is bastardly, but for the moment, normalize to one
    interpolation_value /= 1.6;
    // more bastardlyness
    if(interpolation_value > 1.0)
      interpolation_value = 1.0;

    sum_precipitate_tensors += (_Cijkl_precipitates_MP[_qp])*interpolation_value;
    sum_interpolation_values += interpolation_value;
  }

  //more error checking
  if(sum_interpolation_values > 1.0)
    sum_interpolation_values = 1.0;

  //local elasticity tensor
  sum_precipitate_tensors += _Cijkl_MP[_qp]*(1.0 - sum_interpolation_values);
  _elasticity_tensor[_qp] = sum_precipitate_tensors;

  _Jacobian_mult[_qp] = _elasticity_tensor[_qp];
}

void
LinearSingleCrystalPrecipitateMaterial::computeQpEigenstrain()
{
  Real interpolation_value(0.0);
  Real d_interp_value(0.0);
  Real temp;

  for(unsigned int i=0; i<_n_variants; i++)
  {
    //This is bastardly; later don't forget to normalize the OP to one
    temp = (*_coupled_variables[i])[_qp];///1.6;
    interpolation_value = temp*temp;
    //if (interpolation_value > 1.0)
    // interpolation_value = 1.0;

    d_interp_value = 2.0*temp;

    // Fill in the precipitates' eigenstrains materials property
    (_eigenstrains_MP[_qp])[i] = _eigenstrains_rotated[i]*interpolation_value;
    (_d_eigenstrains_MP[_qp])[i] = _eigenstrains_rotated[i]*d_interp_value;
    (_precipitate_eigenstrain[_qp])[i] = _eigenstrains_rotated[i];
  }
}

void
 LinearSingleCrystalPrecipitateMaterial::computeQpElasticStrain()
 {
   // compute the elastic strain: e_el = e_local - e_misfit

   computeQpStrain();

   // sum up the misfit strains for the orientation variants
   RankTwoTensor sum_precipitate_strains;
   sum_precipitate_strains.zero();

   for(unsigned int i=0; i<_n_variants; i++)
     sum_precipitate_strains += (_eigenstrains_MP[_qp])[i];

   _misfit_strain[_qp] = sum_precipitate_strains;

   _elastic_strain[_qp] = _local_strain[_qp] - _misfit_strain[_qp];
}

void
LinearSingleCrystalPrecipitateMaterial::computeQpElasticStress()
{
  computeQpStress();
}

void
LinearSingleCrystalPrecipitateMaterial::computeQpStrain()
{
  //strain = (grad_disp + grad_disp^T)/2
  RankTwoTensor grad_tensor(_grad_disp_x[_qp],_grad_disp_y[_qp],_grad_disp_z[_qp]);

  _local_strain[_qp] = (grad_tensor + grad_tensor.transpose())/2.0;
}

void
LinearSingleCrystalPrecipitateMaterial::computeQpStress()
{
  // stress = C * e
  _stress[_qp] = _elasticity_tensor[_qp]*_elastic_strain[_qp];
}
