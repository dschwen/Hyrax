/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  6 December 2012
*
*  NucleationMarker is basically a clone of RandomHitMarker in moose_test,
*  thanks go to Derek Gaston.
*
*************************************************************************/

#include "NucleationMarker.h"
#include "NucleationLocationUserObject.h"
#include <ostream>

template<>
InputParameters validParams<NucleationMarker>()
{
  InputParameters params = validParams<Marker>();
    params.addRequiredParam<UserObjectName>("nucleation_userobject", "The name of the UserObject that calculates nucleation event positions");
  return params;
}

NucleationMarker::NucleationMarker(const std::string & name, InputParameters parameters) :
    Marker(name, parameters),
    _nucleation_userobject(getUserObject<NucleationLocationUserObject>("nucleation_userobject"))
{
}

Marker::MarkerValue
NucleationMarker::computeElementMarker()
{
  if(_nucleation_userobject.elementWasHit(_current_elem))
    return REFINE;

  return DONT_MARK;
}
