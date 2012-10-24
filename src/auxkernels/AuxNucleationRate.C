/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*************************************************************************/

#include "AuxNucleationRate.h"

#include <cmath>
#include <ostream>

template<>
InputParameters validParams<AuxNucleationRate>()
{
  InputParameters params = validParams<AuxKernel>();
  //params.addRequiredParam<Real>("Kn1", "First nucleation rate coefficient");
  params.addRequiredParam<Real>("Kn2", "Second nucleation rate coefficient");
  params.addRequiredCoupledVar("coupled_aux_var","coupled auxiliary variable: supersaturation");
  // I'd have line in input file, coupled_aux_var = supersaturation
  params.addRequiredParam<Real>("Z", "Non-equilibrium Zeldovitch factor");
  //params.addRequiredParam<Real>("N", "# of atoms in phase field cell");
  params.addRequiredParam<Real>("Beta_star", "1/characteristic nucleation time");
  params.addRequiredParam<Real>("linear_density", "linear atomic density of matrix");

  return params;
}

AuxNucleationRate::AuxNucleationRate(const std::string & name, InputParameters parameters)
  : AuxKernel(name, parameters),
    _coupled_supersaturation(coupledValue("coupled_aux_var")),
    //_Kn1(getParam<Real>("Kn1")),
    _Kn2(getParam<Real>("Kn2")),
    _Z(getParam<Real>("Z")),
    //_N(getParam<Real>("N")),
    _beta_star(getParam<Real>("Beta_star")),
    _linear_density(getParam<Real>("linear_density"))
{
}

Real
AuxNucleationRate::computeValue()
{
 // j_star = Kn1 * exp(-1*Kn2 / supersaturation)
//adaptive meshing needs to be handled here - need to know the size of the element.

  Real kn1;

  // handling the dimension of the problem here and making sure we get the correct
  // (areal or volume) density
  if (_dim == 2)
  {
    kn1 = _Z*_beta_star*pow(_linear_density, 2);
  }
  else if (_dim == 3)
  {
    kn1 = _Z*_beta_star*pow(_linear_density, 3);
  }
  else
    mooseError("honky, your problem dimesion must be 2 or 3 (AuxNucleationRate");

  // std::cout<<"current elem volume = "<<_current_elem_volume<<std::endl;
  //std::cout<<"current elem pointer = "<<_current_elem<<std::endl;

  // correct the density to the actual element volume to get # of atoms
  kn1 *= _current_elem_volume;

  //return _Kn1*exp(-1.0*_Kn2/_coupled_supersaturation[_qp]);
  return kn1*exp(-1.0*_Kn2/_coupled_supersaturation[_qp]);
}

