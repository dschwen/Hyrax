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
  params.addRequiredParam<Real>("x_center", "Center of the seed in the X direction");
  params.addRequiredParam<Real>("y_center", "Center of the seed in the Y direction");
  params.addParam<Real>("z_center", 0.0, "Center of the seed in the Z direction");

  return params;
}


ACNucleus::ACNucleus(const std::string & name, InputParameters parameters)
    :Kernel(name, parameters),
     _n_value(getParam<Real>("n_value")),
     _radius(getParam<Real>("radius")),
     _int_width(getParam<Real>("int_width")),
     _x_center(getParam<Real>("x_center")),
     _y_center(getParam<Real>("y_center")),
     _z_center(getParam<Real>("z_center"))
{

  std::cout << "in ACnucleus constructor";

}

Real
ACNucleus::computeQpResidual()
{
  Real x = _q_point[_qp](0);
  Real y = _q_point[_qp](1);
  Real z = _q_point[_qp](2);

  /* std::cout<< "x-center " << _x_center;
  std::cout<< "y-center " << _y_center;
  std::cout<< "z-center " << _z_center;
  std::cout<< "x " << x;
  std::cout<< "y " << y;
  std::cout<< "z " << z; */

  //determine if the point is within the radius of the circle
  Real distance = 0.0;
  distance = (x - _x_center)*(x - _x_center)
           + (y - _y_center)*(y - _y_center)
           + (z - _z_center)*(z - _z_center);
  distance = sqrt(distance);

  //if(_u[_qp]>=1.5)
  //  return _u[_qp];
  //else
  //{
    if(distance <=_radius - _int_width/2.0)
      return -1.0* _n_value*_test[_i][_qp];
    else if(distance < _radius + _int_width/2.0)
    {
      Real interface_position = (distance - _radius + _int_width/2.0)/_int_width;
      return -1.0*_test[_i][_qp]*(_n_value - _n_value*(1.0+std::cos(interface_position*libMesh::pi))/2.0);
    }
    else
      return 0.0;
    //}
}
