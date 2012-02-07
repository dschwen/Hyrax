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
#include "SymmTensor.h"
#include "SymmAnisotropicElasticityTensor.h"

// may need some more includes here.

/**
 * LinearSingleCrystalPrecipitateMaterial handles anisotropic, single-crystal material elastic
 * constants.  It's designed to work with SymmSingleCrystalMaterial.  It handles a single
 * crystal of matrix with an arbitrary number of orientation variants of a coherent precipitate.
 * The rotation of orientation variants is currently set up around the C-axis (2D rotation).
 */

template<>
InputParameters validParams<LinearSingleCrystalPrecipitateMaterial>()
{
  InputParameters params = validParams<SolidMechanicsMaterial>();
  params.addRequiredParam<std::vector<Real> >("C_matrix", "Stiffness tensor for matrix");  
  params.addRequiredParam<std::vector<Real> >("C_precipitate", "Stiffness tensor for precipitate");
  params.addRequiredParam<std::vector<Real> >("e_precipitate","Eigenstrain tensor for precipitate");
  params.addRequiredParam<int>("n_variants","# of orientation variants for precipitate in single crystal");
  params.addRequiredCoupledVar("variable_names","Array of coupled variable names");  
  
  return params;
}

template <>
PropertyValue *
MaterialProperty<std::vector<SymmTensor> >::init(int size)
{
  typedef MaterialProperty<std::vector<SymmTensor> > PropType;
  PropType *copy = new PropType;
  copy->_value.resize(size);

  // We don't know the size of the underlying vector at each
  // quadrature point, the user will be responsible for resizing it
  // and filling in the entries...

  // Return the copy we allocated
  return copy;
}

template <>
PropertyValue *
MaterialProperty<std::vector<SymmAnisotropicElasticityTensor> >::init(int size)
{
  typedef MaterialProperty<std::vector<SymmAnisotropicElasticityTensor> > PropType;
  PropType *copy = new PropType;
  copy->_value.resize(size);

  // We don't know the size of the underlying vector at each
  // quadrature point, the user will be responsible for resizing it
  // and filling in the entries...

  // Return the copy we allocated
  return copy;
}


LinearSingleCrystalPrecipitateMaterial::LinearSingleCrystalPrecipitateMaterial(const std::string & name, 
                                                                               InputParameters parameters)
    : SolidMechanicsMaterial(name, parameters),
      _Cijkl_matrix_vector(getParam<std::vector<Real> >("C_matrix")),
      _Cijkl_precipitate_vector(getParam<std::vector<Real> >("C_precipitate")),
      _eigenstrain_vector(getParam<std::vector<Real> >("e_precipitate")),
      _n_variants(getParam<int>("n_variants")),
      _Cijkl_matrix(),
      _Cijkl_precipitates_rotated(),
      _eigenstrains_rotated(),
      _local_strain(declareProperty<SymmTensor >("local_strain")),
      _misfit_strain(declareProperty<SymmTensor >("misfit_strain")),
      _eigenstrains_rotated_MP(declareProperty<std::vector<SymmTensor> >("eigenstrains_rotated_MP")),
      // something about this being a SymmAnisotropicElasticityTensor instead of a SymmElasticityTensor is not ok.
      // Also it seems like the std::vector of tensors is not ok either.
      _Cijkl_matrix_MP(declareProperty<SymmAnisotropicElasticityTensor *>("Cijkl_matrix_MP")),
      _Cijkl_precipitates_rotated_MP(declareProperty<std::vector<SymmAnisotropicElasticityTensor > >("Cijkl_precipitates_rotated_MP"))
 {
   // check to make sure the input file is all set up right
   if(_n_variants != coupledComponents("variable_names"))
     mooseError("Please match the number of orientation variants with coupled order parameters.");

   //size the array to the number of coupled order parameters
   _coupled_variables.resize(_n_variants);

   // loop through the coupled order parameters and couple them in
   for(int i=0; i < _n_variants; i++)
     _coupled_variables[i] = &coupledValue("variable_names", i);

   // // fill in the local tensors from the input vector information
   //
   // going to need to override this or something, or make a note about what to feed in the bool...
   // or add the bool up in the input
   // //_Cijkl_matix.fillFromInputVector(_Cijkl_matrix_vector, /*bool*/);

// go through a loop to fill in these two vectors of SymmTensors
   // for(int i=0; i<_n_variants; i++)
   //{
   //_Cijkl_precipitates_rotated[i].fillFromInputVector(_Cijkl_precipitate_vector, /*bool*/);
   ///// this one requires a new method for SymmTensor
   //_eigenstrains_rotated[i].fillFromInputVector(_eigenstrain_vector, /*bool?*/);
   //}
   //
   //
 
   //resize the vectors to hold the rotated information
    _Cijkl_precipitates_rotated.resize(_n_variants);
    _eigenstrains_rotated.resize(_n_variants);

    // rotate all the things
    for(int i=0; i<_n_variants; i++)
    {
      // in degrees, since SymmAnisotropicElasticityTensor takes degrees
      Real rotation_angle_base = 360.0/Real(_n_variants);
      Real rotation_angle = 0.0;

      // set up temp variables for the precipitate rotations
      SymmAnisotropicElasticityTensor C_tensor(/*& new_Cijkl_precipitate*/);
      SymmTensor e_strain(/*& new_eigenstrain*/);
    
      // do the rotation
//      C_tensor.rotate(rotation_angle);
   //   // C_ijkl_precipitate_rotated[i] = C_tensor;

   //   // e_strain.rotate(rotation_angle);
   //   _eigenstrains_rotated[i] = e_strain;
  
   //   // increment the rotation angle for the next go-round
      rotation_angle = rotation_angle + rotation_angle_base; 
    }
 }

 LinearSingleCrystalPrecipitateMaterial::~LinearSingleCrystalPrecipitateMaterial()
 {
   // delete new_Cijkl_matrix;
   // delete new_Cijkl_precipitate;
   // delete new_eigenstrain;
 }

 void
 LinearSingleCrystalPrecipitateMaterial::computeQpProperties()
 {
   computeQpElasticityTensor();
   computeQpElasticStrain();
   computeQpElasticStress();
 }

 void
 LinearSingleCrystalPrecipitateMaterial::computeQpElasticityTensor()
 {
 /**
  * Not making homogeneous modulus approximation between matrix and precipitate. Currently
  * doing a linear sum of stiffness matrices times the order parameter fraction.  Works only
  * if the OP goes between zero and one!  (Forcing it to max at one, right now).
  **/

  // Fill in the matrix stiffness material property (this WILL work fine)
  //_Cijkl_matrix_MP[_qp] = _Cijkl_matrix;

  // Sum the order parameters and stiffnesses for the precipitates
    SymmElasticityTensor sum_precipitate_tensors(0.0);
    Real sum_order_parameters = 0.0;
  
    for(int i=0; i<_n_variants; i++)
    {
     // Fill in the precipitates' stiffnesses materials property
      //(_Cijkl_precipitates_rotated_MP[_qp])[_i] = _eigenstrains_rotated[i];

      // sum_precipitate_tensors += (*_coupled_variables[i])[_qp](_Cijkl_precipitates_rotated_MP[_qp])[_i];
      sum_order_parameters += (*_coupled_variables[_qp])[i];
    }

    // this is bastardly, but for the moment, truncate the sum to = 1 if it's > 1
    if(sum_order_parameters > 1.0)
      sum_order_parameters = 1.0;
 
    _elasticity_tensor[_qp] = sum_precipitate_tensors /*+ (1.0 - sum_order_parameters)*_Cijkl_matrix_MP[_qp]*/;

   // Jacobian multiplier of stress ... hmm..copying from LinearIsotropicMaterial
    _Jacobian_mult[_qp] = _elasticity_tensor[_qp];
 }

 void
 LinearSingleCrystalPrecipitateMaterial::computeQpElasticStrain()
 {
   // // compute the elastic strain: e_el = e_local - e_misfit
   // /*local strain.  Not adding in the 2x shear strains because the stiffness matrix should do that
   //   with the rotation. */
   _local_strain[_qp] = SymmTensor ( _grad_disp_x[_qp](0),
                                     _grad_disp_y[_qp](1),
                                     _grad_disp_z[_qp](2),
                                     0.5*(_grad_disp_x[_qp](1)+ _grad_disp_y[_qp](0)),
                                     0.5*(_grad_disp_y[_qp](2)+ _grad_disp_z[_qp](1)),
                                     0.5*(_grad_disp_z[_qp](0)+ _grad_disp_x[_qp](2)) );

   // // sum up the misfit strains for the orientation variants
   SymmTensor sum_precipitate_strains(0.0);
   for(int i=0; i<_n_variants; i++)
     // sum_precipitate_strains += (*_eigenstrains_rotated_MP[i])[_qp]*(*_coupled_variables[i])[_qp]*(*_coupled_variables[i])[_qp];
     
     _misfit_strain[_qp] = sum_precipitate_strains;
   
   _elastic_strain[_qp] = _local_strain[_qp] - _misfit_strain[_qp];
 }

 void
 LinearSingleCrystalPrecipitateMaterial::computeQpElasticStress()
 {
   // stress = C * e
   // This is computed here so we can see the field, not do a calculation with it
   _stress[_qp] = _elasticity_tensor[_qp]*_elastic_strain[_qp];
 }
