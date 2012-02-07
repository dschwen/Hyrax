/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  11 January 2012
*
*************************************************************************/

#include "ACTransformElasticDF.h"

/**
 * ACTransformElasticDF handles the elastic energy term for a solid-solid transformation
 * in a phase field model, with coupled conserved and non-conserved parameters.  The
 * variable this kernel is working with is an order parameter.  It can handle an
 * arbitrary number of order parameters.
 */

template<>
InputParameters validParams<ACTransformElasticDF>()
{
  InputParameters params = validParams<ACBulk>();
  params.addRequiredCoupledVar("var_names","Array of coupled variable names minus 1");
  // must go from 0 to n-1 at the moment.
  params.addRequiredParam<int>("n_vars", "# of orientation variants minus 1 for precips in single crystal");
  params.addRequiredParam<int>("OP_number","number of the order parameter for this kernel");
  // must start from 0 at the moment.
  
  return params;
}

ACTransformElasticDF::ACTransformElasticDF(const std::string & name, InputParameters parameters)
    : ACBulk(name, parameters),
      _elasticity_tensor(getMaterialProperty<SymmElasticityTensor>("elasticity_tensor")),
      //_eigenstrains_rotated_MP(getMaterialProperty<std::vector<SymmTensor *> >("eigenstrains_rotated_MP")),
      _local_strain(getMaterialProperty<SymmTensor>("local_strain")),
      _n_vars(getParam<int>("n_vars")),
      _OP_number(getParam<int>("OP_number"))
      {
        // Create a vector of the coupled variables and set = 0 the one that the kernel
        // is operating on 
        if(_n_vars != coupledComponents("var_names"))
          mooseError("Please match the number of orientation variants-1 to coupled OPs -1.");

         _coupled_vars.resize(_n_vars+1);

         for(int i=0; i< _n_vars+1; i++)
           _coupled_vars[i] = &coupledValue("var_names", i);
      }
                         

Real
ACTransformElasticDF::computeDFDOP(PFFunctionType type)
{
// elastic strain = homogeneous strain (BC) + heterogeneous strain - misfit strain
// Follow mostly what is in ACGrGrElasticDF:
  
  switch (type)
  {
  case Residual:

    // hokay, so!
    
    Real elastic_term;
    Real misfit_term;

    elastic_term = calculateLocalTerm();
    misfit_term = calculateMisfitTerm();

    return elastic_term + misfit_term;
    // return 0.0;

  case Jacobian:
    
    return 0.0;
  }
  mooseError("Invalid type passed in");
}
                         
Real
ACTransformElasticDF::calculateLocalTerm()
{
  SymmTensor elastic_a;
  // have to be careful on the multiply operator here - put the scalar at the end of the expression
  elastic_a = (_elasticity_tensor[_qp]*_local_strain[_qp])*2.0;

  SymmTensor elastic_b;
//  elastic_b = (_eigenstrains_rotated_MP[_qp])[_OP_number]*_u[_qp];
                      
  // elastic_term = -2.0 * local_strain * elastic_tensor * eigenstrain_OP * eta_OP
  return elastic_a.doubleContraction(elastic_b);
}

Real
ACTransformElasticDF::calculateMisfitTerm()
{
  Real misfit_term = 0.0;
  SymmTensor misfit_a(0.0);
  SymmTensor misfit_b(0.0);

  //misfit_a = 2.0*(*_eigenstrains_rotated_MP[_qp])[_OP_number]*_u[_qp];
  
  // This will loop over any arbitrary number of order parameters
  for (int i = 0; i < _n_vars+1; i++)
  {
    if (i == _OP_number)
    {
      //misfit_b = (_eigenstrains_rotated_MP[_qp])[i]*_u[_qp];
    }
    else
    {
      // i believe you reference coupled_vars this way
      //misfit_b = (_eigenstrains_rotated_MP[_qp])[_i]*(*_coupled_vars[i])[_qp];
    }
                             
    //misfit_term += 2.0*eta_OP*eigenstrain_OP * elastic_tensor * eigenstrain_i * eta_i
    misfit_term += misfit_a.doubleContraction(misfit_b); 
  }

  return misfit_term;
}
                         
