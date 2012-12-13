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
 */

template<>
InputParameters validParams<LinearSingleCrystalPrecipitateMaterial>()
{
  InputParameters params = validParams<LinearElasticMaterial>();
  // matrix material into C_ijkl
  params.addRequiredParam<std::vector<Real> >("C_precipitate", "Stiffness tensor for precipitate");
  params.addRequiredParam<std::vector<Real> >("e_precipitate","Eigenstrain tensor for precipitate: e11, e22, e33, e23, e13, e12");
  params.addRequiredParam<int>("n_variants","# of orientation variants for precipitate in single crystal");
  params.addRequiredCoupledVar("variable_names","Array of coupled variable names");

  return params;
}

LinearSingleCrystalPrecipitateMaterial::LinearSingleCrystalPrecipitateMaterial(const std::string & name, InputParameters parameters)
    : LinearElasticMaterial(name, parameters),
      //_Cijkl_precipitate_vector(getParam<std::vector<Real> >("C_precipitate")),
      _eigenstrain_vector(getParam<std::vector<Real> >("e_precipitate")),
      _n_variants(getParam<int>("n_variants")),
      //_Cijkl_precipitate(),
      _eigenstrain(),
      //_Cijkl_precipitates_rotated(),
      _eigenstrains_rotated(),
      _local_strain(declareProperty<RankTwoTensor >("local_strain")),
      _misfit_strain(declareProperty<RankTwoTensor >("misfit_strain")),
      _eigenstrains_rotated_MP(declareProperty<std::vector<RankTwoTensor> >("eigenstrains_rotated_MP")),
      _Cijkl_matrix_MP(declareProperty<ElasticityTensorR4>("Cijkl_matrix_MP")),
      //_Cijkl_precipitates_rotated_MP(declareProperty<std::vector<ElasticityTensorR4 > >("Cijkl_precipitates_rotated_MP")),
      //_d_elasticity_tensor(declareProperty<std::vector<ElasticityTensorR4> >("d_elasticity_tensor")),
      _d_eigenstrains_rotated_MP(declareProperty<std::vector<RankTwoTensor> >("d_eigenstrains_rotated_MP"))
{
  // check to make sure the input file is all set up right
  if(_n_variants != coupledComponents("variable_names"))
    mooseError("Please match the number of orientation variants with coupled order parameters (LSXPM).");

  //size the array to the number of coupled order parameters
  _coupled_variables.resize(_n_variants);

  // loop through the coupled order parameters and couple them in
  for(unsigned int i=0; i < _n_variants; i++)
    _coupled_variables[i] = &coupledValue("variable_names", i);

  // using _Cijkl as the matrix material
  // _Cijkl_precipitate.fillFromInputVector(_Cijkl_precipitate_vector, _all_21);
  _eigenstrain.fillFromInputVector(_eigenstrain_vector);

  //resize the vectors to hold the rotated information
  //_Cijkl_precipitates_rotated.resize(_n_variants);
  _eigenstrains_rotated.resize(_n_variants);

  // fill in the first variant without rotation
  //_Cijkl_precipitates_rotated[0] = _Cijkl_precipitate;
  _eigenstrains_rotated[0] = _eigenstrain;

  /*std::cout<<"eigenstraines rotated 0"<<std::endl;
  for(int i=1; i<4; i++)
    for (int j=1; j<4; j++)
    std::cout<<"i="<<i<<" j="<<j<<"val="<<_eigenstrains_rotated[0].getValue(i,j)<<std::endl; */

  // rotate all the things, in radians
   Real rotation_angle_base = 2.0*libMesh::pi/Real(_n_variants);
  // Real rotation_angle_base = 360.0/Real(_n_variants);
  Real rotation_angle = rotation_angle_base;

  for(unsigned int i=1; i<_n_variants; i++)
  {
    // do the rotation
    // _Cijkl_precipitates_rotated[i] = _Cijkl_precipitate.rotate(rotation_angle, 0.0, 0.0);
 // // _Cijkl_precipitates_rotated[i] = _Cijkl_precipitate;

     // e_strain.rotate(rotation_angle);
//    _eigenstrains_rotated[i] = _eigenstrain.rotate(rotation_angle, 0.0, 0.0);
    _eigenstrains_rotated[i] = _eigenstrain.rotateXyPlane(rotation_angle);

    /*  std::cout<<"rotation angle="<<rotation_angle<<std::endl;
    std::cout<<"eigenstraines rotated "<<i<<std::endl;
  for(int k=1; k<4; k++)
    for (int j=1; j<4; j++)
    std::cout<<"i="<<k<<" j="<<j<<"val="<<_eigenstrains_rotated[i].getValue(k,j)<<std::endl; */

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
    //_Cijkl_precipitates_rotated_MP[_qp].resize(_n_variants);
    _eigenstrains_rotated_MP[_qp].resize(_n_variants);
    //_d_elasticity_tensor[_qp].resize(_n_variants);
    _d_eigenstrains_rotated_MP[_qp].resize(_n_variants);
    
    computeQpElasticityTensor();
    computeQpEigenstrain();
    computeQpElasticStrain();
    computeQpElasticStress();
  }
  
 }

 void
 LinearSingleCrystalPrecipitateMaterial::computeQpElasticityTensor()
 {
   //assuming homogeneous modulus between precipitate and matrix
   LinearElasticMaterial::computeQpElasticityTensor();

// ignore all this below for the moment
   /**
    * Not making homogeneous modulus approximation between matrix and precipitate.  Works only
    * if the OP goes between zero and one!  (Forcing it to max at one, right now).
    **/

   // Sum the order parameters and stiffnesses for the precipitates
/*   ElasticityTensorR4 sum_precipitate_tensors;
   ElasticityTensorR4 d_sum_precip_tensors;

   sum_precipitate_tensors.zero();
   d_sum_precip_tensors.zero();

   // interpolation function h
   Real temp;
   Real interpolation_value(0.0), sum_interpolation_values(0.0);
   Real d_interp_value(0.0), d_sum_interp_values(0.0);

   // Fill in the matrix stiffness material property
   _Cijkl_matrix_MP[_qp] = _Cijkl;

   for(unsigned int i=0; i<_n_variants; i++)
    {
      // calculate the value of the interpolation function for an OP and this _qp
      temp = (*_coupled_variables[i])[_qp];
      // this is bastardly, but for the moment, normalize to one
      temp /= 1.6;

      interpolation_value = temp*temp;
      d_interp_value = 2.0*temp;

      // Fill in the precipitates' stiffnesses materials property
      (_Cijkl_precipitates_rotated_MP[_qp])[i] = _Cijkl_precipitates_rotated[i];

      sum_precipitate_tensors += (_Cijkl_precipitates_rotated_MP[_qp])[i]*interpolation_value;
      sum_interpolation_values += interpolation_value;

      d_sum_precip_tensors += (_Cijkl_precipitates_rotated_MP[_qp])[i]*d_interp_value;
      d_sum_interp_values += d_interp_value;

      //derivative of local elasticity tensor
      (_d_elasticity_tensor[_qp])[i] = _Cijkl_matrix_MP[_qp]*d_interp_value*(-1.0)
        + (_Cijkl_precipitates_rotated_MP[_qp])[i]*d_interp_value;
    }

   //local elasticity tensor
   _elasticity_tensor[_qp] = sum_precipitate_tensors
     + _Cijkl_matrix_MP[_qp]*(1.0 - sum_interpolation_values);

   // Jacobian multiplier of stress ... hmm..copying from LinearIsotropicMaterial
   _Jacobian_mult[_qp] = _elasticity_tensor[_qp];
*/
 }

void
LinearSingleCrystalPrecipitateMaterial::computeQpEigenstrain()
{
  Real interpolation_value(0.0);
  Real d_interp_value(0.0);

  for(unsigned int i=0; i<_n_variants; i++)
  {
    // can normalize the OP to one if needed
    interpolation_value = (*_coupled_variables[i])[_qp]*(*_coupled_variables[i])[_qp];
    d_interp_value = 2.0*(*_coupled_variables[i])[_qp];

    // Fill in the precipitates' eigenstrains materials property
    (_eigenstrains_rotated_MP[_qp])[i] = _eigenstrains_rotated[i]*interpolation_value;
    (_d_eigenstrains_rotated_MP[_qp])[i] = _eigenstrains_rotated[i]*d_interp_value;
  }
}

void
 LinearSingleCrystalPrecipitateMaterial::computeQpElasticStrain()
 {
   // compute the elastic strain: e_el = e_local - e_misfit

   //local strain:
   RankTwoTensor grad_tensor(_grad_disp_x[_qp],_grad_disp_y[_qp],_grad_disp_z[_qp]);
   _local_strain[_qp] = (grad_tensor + grad_tensor.transpose())/2.0;

   // // sum up the misfit strains for the orientation variants
   RankTwoTensor sum_precipitate_strains;
   sum_precipitate_strains.zero();

   for(unsigned int i=0; i<_n_variants; i++)
     sum_precipitate_strains += (_eigenstrains_rotated_MP[_qp])[i];

   _misfit_strain[_qp] = sum_precipitate_strains;

   _elastic_strain[_qp] = _local_strain[_qp] - _misfit_strain[_qp];
}

void
LinearSingleCrystalPrecipitateMaterial::computeQpElasticStress()
{
  LinearElasticMaterial::computeQpStress();
}
