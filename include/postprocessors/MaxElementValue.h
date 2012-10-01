/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  1 October 2012
*
*************************************************************************/

#ifndef MAXELEMENTVALUE_H
#define MAXELEMENTVALUE_H

#include "ElementIntegral.h"

class MaxElementValue;

template<>
InputParameters validParams<MaxElementValue>();

class MaxElementValue : public ElementIntegral
{
public:
  MaxElementValue(const std::string & name, InputParameters parameters);

  //virtual void initialize();
  //virtual void execute();
  virtual Real getValue();
  //virtual void threadJoin(const UserObject &y);

protected:

private:
};

#endif //MAXELEMENTVALUE_H
