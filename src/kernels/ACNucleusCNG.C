/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  29 May 2012
*
*************************************************************************/

#include "ACNucleusCNG.h"
#include <cmath>
#include <ostream>

template<>
InputParameters validParams<ACNucleusCNG>()
{
  InputParameters params = validParams<ACNucleus>();

  return params;
}

ACNucleusCNG::ACNucleusCNG(const std::string & name, InputParameters parameters)
    : ACNucleus(name, parameters),
      _nucleation_locations(getMaterialProperty<std::vector<RealVectorValue> >("nucleation_locations")),
      _start_times(getMaterialProperty<std::vector<Real> >("start_times")),
      _end_times(getMaterialProperty<std::vector<Real> >("end_times"))
{
}

Real
ACNucleusCNG::computeQpResidual()
{
  for(int i = 0; i < _nucleation_locations[_qp].size(); i++)
  {
    // determine if the i-th nucleus is active
    if( _t >= _start_times[_qp][i] && _t <= _end_times[_qp][i] )
    {
      // coordinates of current quadrature point
      Real x = _q_point[_qp](0);
      Real y = _q_point[_qp](1);
      Real z = _q_point[_qp](2);

      // determine if distance of current qp from center of a nucleus
      Real distance = 0.0;
      distance = (x - _nucleation_locations[_qp][i](0))*(x - _nucleation_locations[_qp][i](0))
        + (y - _nucleation_locations[_qp][i](1))*(y - _nucleation_locations[_qp][i](1))
        + (z - _nucleation_locations[_qp][i](2))*(z - _nucleation_locations[_qp][i](2));
      distance = sqrt(distance);

      // determine if current qp is in bulk or interface or not at all of nucleus
      if(distance <= _radius - _int_width/2.0)
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
  return 0.0;
}

