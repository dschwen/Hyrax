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
  //params.addRequiredParam<Real>("Kn2", "Second nucleation rate coefficient");

  params.addRequiredCoupledVar("coupled_aux_var","coupled auxiliary variable for free energy change");
  // I'd have line in input file, coupled_aux_var = supersaturation, or chem_elastic, or whatever

  params.addRequiredParam<Real>("gamma", "Surface energy");
  params.addParam<Real>("Kb", 1.3806503e-23, "Boltzmann's constant, make sure units all match");
  params.addParam<Real>("temperature", 473, "Temperature");
  params.addRequiredParam<Real>("scale_factor","factor to scale energy/dimensions by");

  params.addRequiredParam<Real>("Z", "Non-equilibrium Zeldovitch factor");
  //params.addRequiredParam<Real>("N", "# of atoms in phase field cell");
  params.addRequiredParam<Real>("Beta_star", "1/characteristic nucleation time");
  params.addRequiredParam<Real>("linear_density", "linear atomic density of matrix");

  params.addRequiredParam<int>("n_OP_vars", "# of coupled OP variables");
  params.addRequiredCoupledVar("OP_var_names", "Array of coupled OP variable names");

  return params;
}

AuxNucleationRate::AuxNucleationRate(const std::string & name, InputParameters parameters)
  : AuxKernel(name, parameters),
    _coupled_energy(coupledValue("coupled_aux_var")),
    //_Kn1(getParam<Real>("Kn1")),
    // _Kn2(getParam<Real>("Kn2")),
    _Z(getParam<Real>("Z")),
    //_N(getParam<Real>("N")),
    _beta_star(getParam<Real>("Beta_star")),
    _linear_density(getParam<Real>("linear_density")),
    _n_OP_vars(getParam<int>("n_OP_vars")),
    _gamma(getParam<Real>("gamma")),
    _Kb(getParam<Real>("Kb")),
    _temperature(getParam<Real>("temperature")),
    _scale_factor(getParam<Real>("scale_factor"))
{
  if(_n_OP_vars != coupledComponents("OP_var_names"))
    mooseError("Please match the number of orientation variants to coupled OPs (AuxNucleationRate).");

  _coupled_OP_vars.resize(_n_OP_vars);

  for(unsigned int i=0; i< _n_OP_vars; i++)
    _coupled_OP_vars[i] = &coupledValue("OP_var_names", i);
}

Real
AuxNucleationRate::computeValue()
{
  Real kn1;
  Real kn2;
  Real alpha;

  // handling the dimension of the problem here and making sure we get the correct
  // (areal or volume) density
  if (_dim == 2)
  {
    kn1 = _Z*_beta_star*pow(_linear_density, 2);
    alpha = libMesh::pi;
  }
  else if (_dim == 3)
  {
    kn1 = _Z*_beta_star*pow(_linear_density, 3);
    alpha = (16*libMesh::pi)/3;
  }
  else
    mooseError("honky, your problem dimesion must be 2 or 3 (AuxNucleationRate");

  // correct the density to the actual element volume to get # of atoms
  kn1 *= _current_elem_volume;

  kn2 = _scale_factor*alpha*std::pow(_gamma, (int)_dim)/(_Kb*_temperature);

  // We check to see if we're in a particle in AuxNucleationProbability.
  /* for(unsigned int i=0; i<_n_OP_vars; i++)
  {
    if((*_coupled_OP_vars[i])[_qp] > 0.1)
    {
      return 0.0;
    }

    }*/

  return kn1*std::exp(-1.0*kn2/std::pow(_coupled_energy[_qp], (int)_dim-1));
}

