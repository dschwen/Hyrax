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
#include <ostream>

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
  params.addRequiredCoupledVar("OP_var_names","Array of coupled variable names");
params.addRequiredParam<int>("n_OP_vars", "# of orientation variants for precips in single crystal");
  params.addRequiredParam<int>("OP_number","# of the order parameter for this kernel, starting from 1");

  return params;
}

ACTransformElasticDF::ACTransformElasticDF(const std::string & name, InputParameters parameters)
    : ACBulk(name, parameters),
      _elasticity_tensor(getMaterialProperty<SymmElasticityTensor>("elasticity_tensor")),
      _eigenstrains_rotated_MP(getMaterialProperty<std::vector<SymmTensor > >("eigenstrains_rotated_MP")),
      _local_strain(getMaterialProperty<SymmTensor>("local_strain")),
      _n_OP_vars(getParam<int>("n_OP_vars")),
      _OP_number(getParam<int>("OP_number"))
      {
        // Create a vector of the coupled variables and set = 0 the one that the kernel
        // is operating on 
        if(_n_OP_vars != coupledComponents("OP_var_names"))
          mooseError("Please match the number of orientation variants to coupled OPs.");

         _coupled_vars.resize(_n_OP_vars);

         for(int i=0; i< _n_OP_vars; i++)
         {
           if(i == _OP_number-1)
           {
            _coupled_vars[i] = NULL;
           }
           else
           {
           _coupled_vars[i] = &coupledValue("OP_var_names", i);
           }
         }
      }


Real
ACTransformElasticDF::computeDFDOP(PFFunctionType type)
{
// elastic strain = homogeneous strain (BC) + heterogeneous strain - misfit strain
// Follow mostly what is in ACGrGrElasticDF:

  
  Real elastic_term;
  Real misfit_term;
  
  switch (type)
  {
  case Residual:
    
    // hokay, so!
    elastic_term = calculateLocalTerm();
    // std::cout << "elastic term", std::cout << elastic_term, std::cout << std::endl;
    
    misfit_term = calculateMisfitTerm();
    //std::cout << "misfit term", std::cout << misfit_term, std::cout  << std::endl;
    
    return elastic_term + misfit_term;
    
  case Jacobian:
    elastic_term = calculateLocalJacobianTerm();
    misfit_term = calculateMisfitJacobianTerm();

    return 0.0;
    
    //return (elastic_term + misfit_term)*_phi[_j][_qp];
  }
  mooseError("Invalid type passed in");
}

Real
ACTransformElasticDF::calculateLocalTerm()
{
  SymmTensor elastic_a;
  // have to be careful on the multiply operator here - put the scalar at the end of the expression
  elastic_a = (_elasticity_tensor[_qp]*_local_strain[_qp])*(-2.0);

//  std::cout << "elastic_a" << std::cout << elastic_a;

  SymmTensor elastic_b;
  elastic_b = (_eigenstrains_rotated_MP[_qp])[_OP_number]*_u[_qp];
//   std::cout << "elastic_b" << std::cout << elastic_b;
  

  // elastic_term = -2.0 * local_strain * elastic_tensor * eigenstrain_OP * eta_OP
  return elastic_a.doubleContraction(elastic_b);
}

Real
ACTransformElasticDF::calculateMisfitTerm()
{
  Real misfit_term = 0.0;
  SymmTensor misfit_a(0.0);
  SymmTensor misfit_b(0.0);

 misfit_a = (_eigenstrains_rotated_MP[_qp])[_OP_number]*_u[_qp]*2.0;
 misfit_a = _elasticity_tensor[_qp]*misfit_a;
 

  // This will loop over any arbitrary number of order parameters
  for (int i = 0; i < _n_OP_vars; i++)
  {
    if (i == _OP_number-1)
    {
      misfit_b = (_eigenstrains_rotated_MP[_qp])[i]*_u[_qp]*_u[_qp];
    }
    else
    {
      misfit_b = (_eigenstrains_rotated_MP[_qp])[i]*(*_coupled_vars[i])[_qp];
      misfit_b *= (*_coupled_vars[i])[_qp];
    }

    //misfit_term += 2.0*eta_OP*eigenstrain_OP * elastic_tensor * eigenstrain_i * eta_i
    misfit_term += misfit_a.doubleContraction(misfit_b);
  }

  return misfit_term;
}

 Real
 ACTransformElasticDF::calculateLocalJacobianTerm()
 {
   SymmTensor elastic_a;
   elastic_a = (_elasticity_tensor[_qp]*_local_strain[_qp])*(-2.0);

   SymmTensor elastic_b;
   elastic_b = (_eigenstrains_rotated_MP[_qp])[_OP_number];

   return elastic_a.doubleContraction(elastic_b);
 }

 Real
 ACTransformElasticDF::calculateMisfitJacobianTerm()
 {
   Real misfit_term = 0.0;
   SymmTensor misfit_a(0.0);
   SymmTensor misfit_b(0.0);

   misfit_a = (_eigenstrains_rotated_MP[_qp])[_OP_number]*2.0;
   misfit_a = _elasticity_tensor[_qp]*misfit_a;

   for (int i = 0; i < _n_OP_vars; i++)
   {
     if (i == _OP_number-1)
     {
       misfit_b = (_eigenstrains_rotated_MP[_qp])[i]*_u[_qp]*_u[_qp]*3.0;
     }
     else
     {
       misfit_b = (_eigenstrains_rotated_MP[_qp])[i]*(*_coupled_vars[i])[_qp];
       misfit_b *= (*_coupled_vars[i])[_qp];
     }

     misfit_term += misfit_a.doubleContraction(misfit_b);
   }

   return misfit_term;
 }
