/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  11 January 2012
*
*************************************************************************/

#ifndef ACTRANSFORMELASTICDF_H
#define ACTRANSFORMELASTICDF_H

#include "ACBulk.h"

#include "SymmElasticityTensor.h"
#include "SymmTensor.h"

#include <string>

/**
 * ACTransformElasticDF handles the elastic energy term for a solid-solid transformation
 * in a phase field model, with coupled conserved and non-conserved parameters.  The
 * variable this kernel is working with is an order parameter.  It can handle an
 * arbitrary number of order parameters.
 */

// Forward declaration
class ACTransformElasticDF;

template<>
InputParameters validParams<ACTransformElasticDF>();

/**
 * ACTransformElasticDF handles the elastic energy term for a solid-solid transformation
 * in a phase field model, with coupled conserved and non-conserved parameters.
 */
class ACTransformElasticDF : public ACBulk
{
public:

  ACTransformElasticDF(const std::string & name, InputParameters parameters);

protected:

  /**
   * computeDFDOP()
   * @return returns the partial(elastic free energy/order parameter)
   */
  virtual Real computeDFDOP(PFFunctionType type);

  /**
   * calculateLocalTerm()
   * @return returns the portion of the partial derivative containing the local strain term.
   */
  virtual Real calculateLocalTerm();
  
  /**
   * calculateMisfitTerm()
   * @return returns the portion of the partial derivative containing the misfit strain term.
   * This can handle an arbitrary number of order parameters.
   */
virtual Real calculateMisfitTerm();

/**
   * calculateLocalJacobianTerm()
   * @return returns the portion of the partial derivative Jacobian containing the local
   * strain term.
   */
  virtual Real calculateLocalJacobianTerm();
  
  /**
   * calculateMisfitJacobianTerm()
   * @return returns the portion of the partial derivative Jacobian containing the misfit strain term.
   * This can handle an arbitrary number of order parameters.
   */
virtual Real calculateMisfitJacobianTerm();
  
private:

  // system elasticity tensor, varies in space
  MaterialProperty<SymmElasticityTensor> & _elasticity_tensor;
  MaterialProperty<std::vector<SymmTensor > > & _eigenstrains_rotated_MP;
  MaterialProperty<SymmTensor> & _local_strain;

  // number of orientation variants - 1  
  int _n_vars;

  // orientation variant number for this kernel (0 to n-1)
  int _OP_number;

  // Vector of references to the coupled order parameters
  std::vector<VariableValue *> _coupled_vars;
  
};

#endif //ACTRANSFORMELASTICDF_H
