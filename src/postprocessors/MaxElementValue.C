/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  1 October 2012
*
*************************************************************************/

#include "MaxElementValue.h"

template<>
InputParameters validParams<MaxElementValue>()
{
  InputParameters params = validParams<ElementIntegral>();

  return params;
}

MaxElementValue::MaxElementValue(const std::string & name, InputParameters parameters) :
    ElementIntegral(name, parameters)
{
}

Real
MaxElementValue::getValue()
{
  // was returning same value as ElementIntegral.  Going to scrap this
  // PP for the moment.

  //gatherMax(_integral_value);
  //return _integral_value;
  return 0.0;
}
