/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  1 November 2013
*
*************************************************************************/

#include "Heat.h"

template<>
InputParameters validParams<Heat>()
{
  InputParameters params = validParams<Diffusion>();

  //any coupled variables to add?

  return params;
}


Heat::Heat(const std::string & name, InputParameters parameters) :
    Diffusion(name, parameters),
    _diffusivity(getMaterialProperty<Real>("thermal_diffusivity")),
    _dDiffusivity_dT(getMaterialProperty<Real>("dThermal_diffusivity_dT"))
    //any coupled variables to add?
{
}

Real
Heat::computeQpResidual()
{
  return _grad_test[_i][_qp]*_diffusivity[_qp]*_grad_u[_qp];
}

Real
Heat::computeQpJacobian()
{
  //this might not be perfect, but uh.  Let's go with it.
  return _grad_phi[_j][_qp]*(_dDiffusivity_dT[_qp]*_grad_u[_qp]
                              + _diffusivity[_qp]*_grad_test[_i][_qp]);
}
