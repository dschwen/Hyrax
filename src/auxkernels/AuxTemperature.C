/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  11 October 2013
*
*************************************************************************/

#include "AuxTemperature.h"

template<>
InputParameters validParams<AuxTemperature>()
{
  InputParameters params = validParams<AuxKernel>();

  params.addParam<Real>("temp_in_K", 600, "Simulation temperature in K");

  return params;
}

AuxTemperature::AuxTemperature(const std::string & name, InputParameters parameters)
    : AuxKernel(name, parameters),
      _T(getParam<Real>("temp_in_K"))
{
  std::cout<<"AuxTemperature constructor"<<std::endl;
}

Real
AuxTemperature::computeValue()
{
  return _T;
}

