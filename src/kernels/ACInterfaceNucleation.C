/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  12 June 2012
*
*************************************************************************/

#include "ACInterfaceNucleation.h"

template<>
InputParameters validParams<ACInterfaceNucleation>()
{
  InputParameters params = validParams<ACInterface>();
  params.addParam<Real>("start_time", 0.0, "when to start operating this kernel");
  params.addParam<Real>("end_time", 0.1, "when to stop operating this kernel");
  params.addRequiredParam<Real>("radius", "the radius of the seed circle");
//  params.addParam<Real>("int_width",0.0, "Interface width of the seed circle, defaults to 0");
  params.addParam<RealVectorValue>("center_point", RealVectorValue(0,0,0), "Center of the seed in the X, y, and z directions");

  return params;
}

ACInterfaceNucleation::ACInterfaceNucleation(const std::string & name, InputParameters parameters)
    : ACInterface(name, parameters),
      _start_time(getParam<Real>("start_time")),
      _end_time(getParam<Real>("end_time")),
      _radius(getParam<Real>("radius")),
//      _int_width(getParam<Real>("int_width")),
      _nucleation_center(getParam<RealVectorValue>("center_point"))//,
{
}

RealGradient
ACInterfaceNucleation::precomputeQpResidual()
{
  // if during nucleation event, see what to do
  // else do normal
  if(_t >= _start_time && _t < _end_time)
  {
    // if outside of the seed, do normal
    // if inside of the seed, return 0.0
     Point current_point;
     current_point = _q_point[_qp];

     Real distance;
     distance = (current_point - _nucleation_center).size();

     // determine if current qp is in nucleus or not
     // if(distance > _radius + _int_width/2.0)
     if(distance > _radius)
       return ACInterface::precomputeQpResidual();
     else
       return 0.0;
  }
  else
    return ACInterface::precomputeQpResidual();
}

RealGradient
ACInterfaceNucleation::precomputeQpJacobian()
{
  // if during nucleation event, see what to do
  // else do normal
  if(_t >= _start_time && _t < _end_time)
  {
    // if outside of the seed, do normal
    // if inside of the seed, return 0.0
    Point current_point;
    current_point = _q_point[_qp];

    Real distance;
    distance = (current_point - _nucleation_center).size();

    // determine if current qp is in nucleus or not
    //if(distance > _radius + _int_width/2.0)
    if(distance > _radius)
      return ACInterface::precomputeQpJacobian();
    else
      return 0.0;
  }
  else
    return ACInterface::precomputeQpJacobian();

  return 0.0;
}


