/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  17 May 2012
*
*************************************************************************/

#include "ACNucleus.h"
#include <cmath>
#include <ostream>

template<>
InputParameters validParams<ACNucleus>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<Real>("n_value", "Order parameter value within the seed");
  params.addRequiredParam<Real>("radius", "the radius of the seed circle");
  params.addParam<Real>("int_width",0.0, "Interface width of the seed circle, defaults to 0");
  params.addParam<Real>("x_center", 0.0, "Center of the seed in the X direction");
  params.addParam<Real>("y_center", 0.0, "Center of the seed in the Y direction");
  params.addParam<Real>("z_center", 0.0, "Center of the seed in the Z direction");
  //add something in here about when to operate this
  params.addParam<Real>("start_time", 0.0, "when to start operating this kernel");
  params.addParam<Real>("end_time", 0.1, "when to stop operating this kernel");

  return params;
}


ACNucleus::ACNucleus(const std::string & name, InputParameters parameters)
    :Kernel(name, parameters),
     _radius(getParam<Real>("radius")),
     _int_width(getParam<Real>("int_width")),
     _n_value(getParam<Real>("n_value")),
     _x_center(getParam<Real>("x_center")),
     _y_center(getParam<Real>("y_center")),
     _z_center(getParam<Real>("z_center")),
     _start_time(getParam<Real>("start_time")),
     _end_time(getParam<Real>("end_time"))
{
}

Real
ACNucleus::computeQpResidual()
{
  if(_t >= _start_time && _t < _end_time)
  {
    Real x = _q_point[_qp](0);
    Real y = _q_point[_qp](1);
    Real z = _q_point[_qp](2);

    //determine if the point is within the radius of the circle
    Real distance = 0.0;
    distance = (x - _x_center)*(x - _x_center)
      + (y - _y_center)*(y - _y_center)
      + (z - _z_center)*(z - _z_center);
    distance = sqrt(distance);

    if(distance <=_radius - _int_width/2.0)
      return -1.0* _n_value*_test[_i][_qp];
    else if(distance < _radius + _int_width/2.0)
    {
      Real interface_position = (distance - _radius + _int_width/2.0)/_int_width;
      return -1.0*_test[_i][_qp]*(_n_value - _n_value*(1.0+std::cos(interface_position*libMesh::pi))/2.0);
    }
    else
      return 0.0;
  }
  else
    return 0.0;
}
