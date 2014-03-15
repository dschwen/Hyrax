/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*************************************************************************/

#include "AuxFullNucleationRate.h"

#include <cmath>

template<>
InputParameters validParams<AuxFullNucleationRate>()
{
  InputParameters params = validParams<AuxKernel>();

  params.addRequiredCoupledVar("coupled_bulk_energy_change","coupled auxiliary variable for free energy change");
  // I'd have line in input file, coupled_aux_var = supersaturation, or chem_elastic, or whatever

  params.addRequiredParam<Real>("gamma", "Surface energy");
  params.addParam<Real>("Kb", 1.3806503e-23, "Boltzmann's constant, make sure units all match");
  //params.addParam<Real>("temperature", 473, "Temperature");
  params.addParam<Real>("scale_factor", 1, "factor to scale nucleation rate by");

  //params.addRequiredParam<Real>("Z", "Non-equilibrium Zeldovitch factor");
  //params.addRequiredParam<Real>("Beta_star", "1/characteristic nucleation time");
  params.addRequiredParam<Real>("linear_density", "linear atomic density of matrix");

  params.addRequiredCoupledVar("T", "temperature");
  params.addRequiredCoupledVar("X", "atomic fraction of solute");

  params.addRequiredParam<Real>("jump_distance", "atomic distance");
//  params.addRequiredParam<int>("n_OP_variables", "# of coupled OP variables, >=1");
  // params.addRequiredCoupledVar("OP_variable_names", "Array of coupled OP variable names");


  return params;
}

AuxFullNucleationRate::AuxFullNucleationRate(const std::string & name, InputParameters parameters)
  : AuxKernel(name, parameters),
    _mesh_dimension(_mesh.dimension()),
    _coupled_energy(coupledValue("coupled_bulk_energy_change")),
    _Z(),
    _beta_star(),
    _linear_density(getParam<Real>("linear_density")),
    _gamma(getParam<Real>("gamma")),
    _Kb(getParam<Real>("Kb")),
    _scale_factor(getParam<Real>("scale_factor")),
    //  _n_OP_variables(getParam<int>("n_OP_variables")),
    _T(coupledValue("T")),
    _X(coupledValue("X")),
    _D(getMaterialProperty<Real>("D_alpha")),
    _jump_distance(getParam<Real>("jump_distance")),
    _Omega(getMaterialProperty<Real>("molar_volume"))
{
/*  // Create a vector of the coupled OP variables and gradients
  if(_n_OP_variables != coupledComponents("OP_variable_names"))
    mooseError("Please match the # of orientation variants to coupled OPs (CHCoupledCalphad)");

  _OP.resize(_n_OP_variables);

  for(unsigned int i=0; i< _n_OP_variables; i++)
    _OP[i] = &coupledValue("OP_variable_names", i);
*/
}

Real
AuxFullNucleationRate::computeValue()
{
  computeCriticalRadius();
  computeCriticalEnergy();
  computeZeldovichFactor();
  computeCriticalFrequency();
  computeNumAtoms();

//  std::cout<<"ZNBetastar = "<< _Z*_N*_beta_star <<std::endl;
//  std::cout<<"(-1*G_star)/(Kb*T) = "<<(-1*_G_star)/ (_Kb*_T[_qp]) <<std::endl;
//  std::cout<<"exp(stuff) = "<<std::exp((-1*_G_star)/ (_Kb*_T[_qp]))<<std::endl;

//  std::cout<<"J* scaled ="<<_scale_factor*( _Z*_N*_beta_star*std::exp( (-1*_G_star)/ (_Kb*_T[_qp]) ))<<std::endl;

  return  _scale_factor*( _Z*_N*_beta_star*std::exp( (-1*_G_star)/ (_Kb*_T[_qp]) ));
}

void
AuxFullNucleationRate::computeCriticalRadius()
{
//  if (_mesh_dimension == 2)
//  _r_star = _gamma/(-4*_coupled_energy[_qp]);

  //else if (_mesh_dimension == 3)
  //removed the negative sign
    _r_star = 2*_gamma/_coupled_energy[_qp];

    //else
    //mooseError("honky, your problem dimesion must be 2 or 3 (AuxFullNucleationRate");

//  std::cout<<"r* = "<<_r_star<<std::endl;
}

void
AuxFullNucleationRate::computeCriticalEnergy()
{
  // if (_mesh_dimension == 2)
  //_G_star = (4*libMesh::pi*_r_star*_r_star*_coupled_energy[_qp])
  //  + 2*libMesh::pi*_r_star*_gamma;

  //else if (_mesh_dimension == 3)
  // _G_star = (4*libMesh::pi*_r_star*_r_star*_r_star*_coupled_energy[_qp])/3
  //  + 4*libMesh::pi*_r_star*_r_star*_gamma;

    //else
    //mooseError("honky, your problem dimesion must be 2 or 3 (AuxFullNucleationRate");
  //unsigned int d = (int) _mesh_dimension;

 Real alpha;

 // if (_mesh_dimension == 2)
 // alpha = libMesh::pi;

 // else if (_mesh_dimension == 3)
    alpha = 16*libMesh::pi/3;

    // else
    //mooseError("honky, your problem dimesion must be 2 or 3 (AuxFullNucleationRate");


  _G_star = alpha*std::pow(_gamma, 3)/std::pow(_coupled_energy[_qp], 2);

  // std::cout<<"G* in Joules = "<<_G_star<<std::endl;

  //_G_star = (4*libMesh::pi*_r_star*_r_star*_coupled_energy[_qp])
  //  + 2*libMesh::pi*_r_star*_gamma;

  //else if (_mesh_dimension == 3)
  // _G_star = (4*libMesh::pi*_r_star*_r_star*_r_star*_coupled_energy[_qp])/3
  //  + 4*libMesh::pi*_r_star*_r_star*_gamma;

    //else
    //mooseError("honky, your problem dimesion must be 2 or 3 (AuxFullNucleationRate");

  // Real critical_volume = 4*libMesh::pi*_r_star*_r_star*_r_star/3;

  //rstar seems kind of large, so I'm going to chop it in half, hence /8)
  //Real critical_volume = 4*libMesh::pi*_r_star*_r_star*_r_star/3/8;

  //here I'm assuming that the molar volume is constant across phases.. meh
  //Real num_atoms = critical_volume*6.02214E23/_Omega[_qp];

  // _G_star /= num_atoms;


//  std::cout<<"critical volume = "<<critical_volume<<std::endl;
//  std::cout<<"num_atoms = "<<num_atoms<<std::endl;
//  std::cout<<"G* = "<<_G_star<<std::endl;
}

void
AuxFullNucleationRate::computeZeldovichFactor()
{
  Real critical_volume = 4*libMesh::pi*_r_star*_r_star*_r_star/3;

  //here I'm assuming that the molar volume is constant across phases.. meh
  Real Nc = critical_volume*6.02214E23/_Omega[_qp];

  _Z =std::sqrt( _G_star/( 3*libMesh::pi*Nc*Nc*_Kb*_T[_qp] ));

  //std::cout<<"Z = "<<_Z<<std::endl;
}

void
AuxFullNucleationRate::computeCriticalFrequency()
{
  Real Zc = 4*libMesh::pi*_r_star*_r_star*_linear_density*_linear_density;

  // Real Zc = 4*libMesh::pi*_r_star*_r_star*_linear_density*_linear_density/4;
// std::cout<<"Zc*X = "<<Zc*_X[_qp]<<std::endl;

  _beta_star = Zc*_X[_qp]*_D[_qp]/ ( std::pow(_jump_distance,2)) ;

//  std::cout<<"beta* = "<<_beta_star<<std::endl;
}

void
AuxFullNucleationRate::computeNumAtoms()
{
  // if (_mesh_dimension == 2)
  // _N = std::pow(_linear_density, 2);

  //else if (_mesh_dimension == 3)
  //  _N = std::pow(_linear_density, 3);

    //else
    //mooseError("honky, your problem dimesion must be 2 or 3 (AuxNucleationRate");

  // correct the density to the actual element volume to get # of atoms
    //THIS IS CURRENTLY WRONG

    Real atomic_volume = _Omega[_qp]/6.02214E23;

    //this is definitely kind of f'd but it's a good start
    Real vol = std::sqrt(_current_elem_volume);
    vol = vol*vol*vol*1E-9*1E-9*1E-9;

    //  std::cout<<"atomic volume = "<<atomic_volume<<std::endl;
    // std::cout<<"volume = "<<vol<<std::endl;

    _N = vol/atomic_volume;

    // std::cout<<"N = "<<_N<<std::endl;
}



