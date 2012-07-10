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
//#include "SymmTensor.h"
//#include "SymmAnisotropicElasticityTensor.h"

#include <ostream>

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
//  params.addRequiredParam<std::vector<Real> >("C_matrix", "Stiffness tensor for matrix: C11, C12, C13, C22, C23, C33, C44, C55, C66 (for 9 inputs)");
  params.addRequiredParam<std::vector<Real> >("C_precipitate", "Stiffness tensor for precipitate");
  params.addRequiredParam<std::vector<Real> >("e_precipitate","Eigenstrain tensor for precipitate: e11, e22, e33, e23, e13, e12");
  params.addRequiredParam<int>("n_variants","# of orientation variants for precipitate in single crystal");
  params.addRequiredCoupledVar("variable_names","Array of coupled variable names");
  // params.addParam<bool>("all_21", false,"Boolean to indicate if all 21 elastic constants are specified in the input (true) or only nine (false)" );

  return params;
}

LinearSingleCrystalPrecipitateMaterial::LinearSingleCrystalPrecipitateMaterial(const std::string & name, InputParameters parameters)
    : LinearElasticMaterial(name, parameters),
      //  _Cijkl_matrix_vector(getParam<std::vector<Real> >("C_matrix")),
      _Cijkl_precipitate_vector(getParam<std::vector<Real> >("C_precipitate")),
      _eigenstrain_vector(getParam<std::vector<Real> >("e_precipitate")),
      _n_variants(getParam<int>("n_variants")),
      //_all_21(getParam<bool>("all_21")),
      // _Cijkl_matrix(),
      _Cijkl_precipitate(),
      _eigenstrain(),
      _Cijkl_precipitates_rotated(),
      _eigenstrains_rotated(),
      _local_strain(declareProperty<RankTwoTensor >("local_strain")),
      _misfit_strain(declareProperty<RankTwoTensor >("misfit_strain")),
      _eigenstrains_rotated_MP(declareProperty<std::vector<RankTwoTensor> >("eigenstrains_rotated_MP")),
      _Cijkl_matrix_MP(declareProperty<RankFourTensor>("Cijkl_matrix_MP")),
      _Cijkl_precipitates_rotated_MP(declareProperty<std::vector<RankFourTensor > >("Cijkl_precipitates_rotated_MP")),
      _d_elasticity_tensor(declareProperty<std::vector<RankFourTensor> >("d_elasticity_tensor")),
      _d_eigenstrains_rotated_MP(declareProperty<std::vector<RankTwoTensor> >("d_eigenstrains_rotated_MP"))
{
  // check to make sure the input file is all set up right
  if(_n_variants != coupledComponents("variable_names"))
    mooseError("Please match the number of orientation variants with coupled order parameters.");

  //size the array to the number of coupled order parameters
  _coupled_variables.resize(_n_variants);

  // loop through the coupled order parameters and couple them in
  for(unsigned int i=0; i < _n_variants; i++)
    _coupled_variables[i] = &coupledValue("variable_names", i);

  // fill in the local tensors from the input vector information
//  _Cijkl_matrix.fillFromInputVector(_Cijkl_matrix_vector, _all_21);

  // using _Cijkl as the matrix material
  _Cijkl_precipitate.fillFromInputVector(_Cijkl_precipitate_vector, _all_21);
  _eigenstrain.fillFromInputVector(_eigenstrain_vector);

  //resize the vectors to hold the rotated information
  _Cijkl_precipitates_rotated.resize(_n_variants);
  _eigenstrains_rotated.resize(_n_variants);

  // fill in the first variant without rotation
  _Cijkl_precipitates_rotated[0] = _Cijkl_precipitate;
  _eigenstrains_rotated[0] = _eigenstrain;

 // rotate all the things, in radians
  Real rotation_angle_base = 2.0*libMesh::pi/Real(_n_variants);
  Real rotation_angle = rotation_angle_base;
//std::cout << "rotation angle base " << rotation_angle_base << std::endl;
//  std::cout << "rotation angle " << rotation_angle << std::endl;

  for(unsigned int i=1; i<_n_variants; i++)
  {
    // do the rotation
    // C_tensor.rotate(rotation_angle, 0.0, 0.0);
     _Cijkl_precipitates_rotated[i] = _Cijkl_precipitate.rotate(rotation_angle, 0.0, 0.0);

    //  _Cijkl_precipitates_rotated[i] = _Cijkl_precipitate;

     // e_strain.rotate(rotation_angle);
    _eigenstrains_rotated[i] = _eigenstrain.rotate(rotation_angle, 0.0, 0.0);
    // increment the rotation angle for the next go-round
    rotation_angle = rotation_angle + rotation_angle_base;
    // std::cout << "rotation angle " << rotation_angle << std::endl;

  }
}

void
 LinearSingleCrystalPrecipitateMaterial::computeQpProperties()
 {
   // resize all the material properties vectors.  Don't forget this.
   _Cijkl_precipitates_rotated_MP[_qp].resize(_n_variants);
   _eigenstrains_rotated_MP[_qp].resize(_n_variants);
   _d_elasticity_tensor[_qp].resize(_n_variants);
   _d_eigenstrains_rotated_MP[_qp].resize(_n_variants);

   computeQpElasticityTensor();
   computeQpEigenstrain();
   computeQpElasticStrain();
   computeQpElasticStress();
 }

 void
 LinearSingleCrystalPrecipitateMaterial::computeQpElasticityTensor()
 {
 /**
  * Not making homogeneous modulus approximation between matrix and precipitate.  Works only
  * if the OP goes between zero and one!  (Forcing it to max at one, right now).
  **/

   // Sum the order parameters and stiffnesses for the precipitates
   RankFourTensor sum_precipitate_tensors;
   RankFourTensor d_sum_precip_tensors;

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
      // this is bastardly, but for the moment, truncate to 1 if > 1
      if (temp > 1.0)
        temp = 1.0;

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

 }

void
LinearSingleCrystalPrecipitateMaterial::computeQpEigenstrain()
{
  Real temp(0.0);
  Real interpolation_value(0.0);
  // Real sum_interpolation_values(0.0);
  Real d_interp_value(0.0);
  //Real d_sum_interp_values(0.0);

  for(unsigned int i=0; i<_n_variants; i++)
  {
    // calculate the value of the interpolation function for an OP and this _qp
    temp = (*_coupled_variables[i])[_qp];

    // this is bastardly, but for the moment, truncate to 1 if > 1
    if (temp > 1.0)
      temp = 1.0;

    interpolation_value = temp*temp;
    d_interp_value = 2.0*temp;

    // Fill in the precipitates' eigenstrains materials property
    (_eigenstrains_rotated_MP[_qp])[i] = _eigenstrains_rotated[i]*interpolation_value;
    (_d_eigenstrains_rotated_MP[_qp])[i] = _eigenstrains_rotated[i]*d_interp_value;
  }
}

void
 LinearSingleCrystalPrecipitateMaterial::computeQpElasticStrain()
 {
   // // compute the elastic strain: e_el = e_local - e_misfit

   //local strain:
   _local_strain[_qp].setValue(0.5*(_grad_disp_x[_qp](0) +
                                    _grad_disp_x[_qp](0)), 1, 1);
   _local_strain[_qp].setValue(0.5*(_grad_disp_y[_qp](1) +
                                    _grad_disp_y[_qp](1)), 2, 2);
   _local_strain[_qp].setValue(0.5*(_grad_disp_z[_qp](2) +
                                    _grad_disp_z[_qp](2)), 3, 3);
   _local_strain[_qp].setValue(0.5*(_grad_disp_x[_qp](1) +
                                    _grad_disp_y[_qp](0)), 1, 2);
   _local_strain[_qp].setValue(0.5*(_grad_disp_x[_qp](1) +
                                    _grad_disp_y[_qp](0)), 2, 1);
   _local_strain[_qp].setValue(0.5*(_grad_disp_x[_qp](2) +
                                    _grad_disp_z[_qp](0)), 1, 3);
   _local_strain[_qp].setValue(0.5*(_grad_disp_x[_qp](2) +
                                    _grad_disp_z[_qp](0)), 3, 1);
   _local_strain[_qp].setValue(0.5*(_grad_disp_y[_qp](2) +
                                    _grad_disp_z[_qp](1)), 2, 3);
   _local_strain[_qp].setValue(0.5*(_grad_disp_y[_qp](2) +
                                    _grad_disp_z[_qp](1)), 3, 2);
   //std::cout << _local_strain[_qp], std::cout << std::endl;

   // // sum up the misfit strains for the orientation variants
   RankTwoTensor sum_precipitate_strains;
   sum_precipitate_strains.zero();

   for(unsigned int i=0; i<_n_variants; i++)
     sum_precipitate_strains += (_eigenstrains_rotated_MP[_qp])[i];

   _misfit_strain[_qp] = sum_precipitate_strains;
   //  std::cout << _misfit_strain[_qp], std::cout << std::endl;

   _elastic_strain[_qp] = _local_strain[_qp] - _misfit_strain[_qp];
   // std::cout << _elastic_strain[_qp], std::cout << std::endl;
 }

 void
 LinearSingleCrystalPrecipitateMaterial::computeQpElasticStress()
 {
   // stress = C * e
   _stress[_qp] = _elasticity_tensor[_qp]*_elastic_strain[_qp];
   //std::cout << _stress[_qp], std::cout << std::endl;
 }
