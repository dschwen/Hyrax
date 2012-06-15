/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  14 June 2012
*
*************************************************************************/

#include "TimeDerivativeNucleation.h"

template<>
InputParameters validParams<TimeDerivativeNucleation>()
{
  InputParameters params = validParams<TimeDerivative>();
  params.addParam<Real>("start_time", 0.0, "when to start operating this kernel");
  params.addParam<Real>("end_time", 0.1, "when to stop operating this kernel");
  params.addRequiredParam<Real>("radius", "the radius of the seed circle");
  params.addParam<RealVectorValue>("center_point", RealVectorValue(0,0,0), "Center of the seed in the X, y, and z directions");

  return params;
}

TimeDerivativeNucleation::TimeDerivativeNucleation(const std::string & name,
                                                   InputParameters parameters)
    : TimeDerivative(name, parameters),
      _start_time(getParam<Real>("start_time")),
      _end_time(getParam<Real>("end_time")),
      _radius(getParam<Real>("radius")),
      _nucleation_center(getParam<RealVectorValue>("center_point"))
{
}

Real
TimeDerivativeNucleation::computeQpResidual()
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

    // Currently, sharp interface only is supported
    if(distance <=_radius)
      return 0.0;
    else
      return TimeDerivative::computeQpResidual();
  }
  else
    return TimeDerivative::computeQpResidual();
}

Real
TimeDerivativeNucleation::computeQpJacobian()
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

    // Currently, sharp interface only is supported
    if(distance <=_radius)
      return 0.0;
    else
      return TimeDerivative::computeQpJacobian();
  }
  else
    return TimeDerivative::computeQpJacobian();
}
