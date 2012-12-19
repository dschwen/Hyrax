/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  19 December 2012
*
*  ErrorFractionMaxHMarker is a clone of ErrorFractionMarker with the
*  added bonus of setting a maximum refinement level.
*
*************************************************************************/

#include "ErrorFractionMaxHMarker.h"

template<>
InputParameters validParams<ErrorFractionMaxHMarker>()
{
  InputParameters params = validParams<ErrorFractionMarker>();

  params.addRequiredParam<unsigned int>("max_h_level", "Maximum h-adapt level for the mesh");
//  params.addParam<Real>("coarsen", 0, "Elements within this percentage of the min error will be coarsened.  Must be between 0 and 1!");
  // params.addParam<Real>("refine", 0, "Elements within this percentage of the max error will be refined.  Must be between 0 and 1!");
  return params;
}


ErrorFractionMaxHMarker::ErrorFractionMaxHMarker(const std::string & name, InputParameters parameters) :
    ErrorFractionMarker(name, parameters),
    //_coarsen(parameters.get<Real>("coarsen")),
    //_refine(parameters.get<Real>("refine"))
    _max_h_level(getParam<unsigned int>("max_h_level"))
{
  //mooseAssert(_coarsen <= 1, "coarsen amount in " + _name + " not less than 1!");
  //mooseAssert(_coarsen >= 0, "coarsen amount in " + _name + " not greater than 0!");
  //mooseAssert(_refine <= 1, "refine amount in " + _name + " not less than 1!");
  //mooseAssert(_refine >= 0, "refine amount in " + _name + " not greater than 0!");
}

/*void
ErrorFractionMarker::markerSetup()
{
  _min = std::numeric_limits<Real>::max();
  _max = 0;

  // First find the max and min error
  for(unsigned int i=0; i<_error_vector.size(); i++)
  {
    _min = std::min(_min, static_cast<Real>(_error_vector[i]));
    _max = std::max(_max, static_cast<Real>(_error_vector[i]));
  }

  _delta = _max-_min;
  _refine_cutoff = (1.0-_refine)*_max;
  _coarsen_cutoff = _coarsen*_delta + _min;
  }*/

Marker::MarkerValue
ErrorFractionMaxHMarker::computeElementMarker()
{
  Real error = _error_vector[_current_elem->id()];

    if((error > _refine_cutoff) && (_current_elem->level() < _max_h_level))
      return REFINE;
    else if(error < _coarsen_cutoff)
      return COARSEN;

  return DO_NOTHING;
}

