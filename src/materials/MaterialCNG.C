/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  30 May 2012
*
*************************************************************************/

#include "MaterialCNG.h"
#include <ostream>

template<>
InputParameters validParams<MaterialCNG>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredCoupledVar("nucleation_field_name", "Same as the variable used for the AuxNucleation auxkernel");
  params.addRequiredParam<Real>("dwell_time", "How long nucleation kernel should be applied to an event");

  return params;
}

MaterialCNG::MaterialCNG(const std::string & name, InputParameters parameters)
    : Material(name, parameters),
      _nucleation_field(coupledValue("nucleation_field_name")),
      _nucleation_locations(declareProperty<std::vector<RealVectorValue> >("nucleation_locations")),
      _start_times(declareProperty<std::vector<Real> >("start_times")),
      _end_times(declareProperty<std::vector<Real> >("end_times")),
      _dwell_time(getParam<Real>("dwell_time"))
{
}

void
MaterialCNG::computeProperties()
{
// The trick here is to resize the vector that holds the nucleation events locations each
// time there is a new event and tack it on to the end.

  int lsize(_locations.size());

  // loop over all quadrature points to get list of nucleation locations
  for(unsigned int qp=0; qp<_qrule->n_points(); ++qp)
  {
    std::cout << "nucleation field[" << qp << "] = " << _nucleation_field[qp]  << std::endl;

    if (_nucleation_field[qp] > 1.0)
   {
      std::cout << "in first for loop if 1" << std::endl;

      // resize the locations vector
      _locations.resize(++lsize);

      // fill in with the point location of the current qp
      _locations[lsize-1](0) = _q_point[qp](0);
      _locations[lsize-1](1) = _q_point[qp](1);
      _locations[lsize-1](2) = _q_point[qp](2);

      // resize the time vectors
      _start.resize(lsize);
      _end.resize(lsize);

      // fill in the time vectors with the start and end times for the new point
      _start[lsize-1] = _t;
      _end[lsize-1] = _t + _dwell_time;
   }
  }

  // loop over all quadrature points again to fill in the materials properties
  for(unsigned int qp=0; qp<_qrule->n_points(); qp++)
  {
    // resize all the materials properties vectors
    _nucleation_locations[qp].resize(lsize);
    _start_times[qp].resize(lsize);
    _end_times[qp].resize(lsize);

    // fill in the vectors
    _nucleation_locations[qp] = _locations;
    _start_times[qp] = _start;
    _end_times[qp] = _end;
  }
}
