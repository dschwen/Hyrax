/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  12 June 2012
*
*************************************************************************/

#include "ACBulkNucleation.h"

template<>
InputParameters validParams<ACBulkNucleation>()
{
  InputParameters params = validParams<ACBulkCoupled>();
  params.addParam<Real>("start_time", 0.0, "when to start operating this kernel");
  params.addParam<Real>("end_time", 0.1, "when to stop operating this kernel");
  params.addRequiredParam<Real>("radius", "the radius of the seed circle");
//  params.addParam<Real>("int_width",0.0, "Interface width of the seed circle, defaults to 0");
  params.addParam<RealVectorValue>("center_point", RealVectorValue(0,0,0), "Center of the seed in the X, y, and z directions");

  return params;
}

ACBulkNucleation::ACBulkNucleation(const std::string & name, InputParameters parameters)
    : ACBulkCoupled(name, parameters),
      _start_time(getParam<Real>("start_time")),
      _end_time(getParam<Real>("end_time")),
      _radius(getParam<Real>("radius")),
//      _int_width(getParam<Real>("int_width")),
      _nucleation_center(getParam<RealVectorValue>("center_point"))//,
{
}

Real
ACBulkNucleation::computeDFDOP(PFFunctionType type)
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
      return ACBulkCoupled::computeDFDOP(type);
    else
      return 0.0;
  }
  else
    return ACBulkCoupled::computeDFDOP(type);
}

