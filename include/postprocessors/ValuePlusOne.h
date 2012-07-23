#ifndef VALUEPLUSONE_H
#define VALUEPLUSONE_H

#include "ChangeVariableData.h"

// forward declaration
class ValuePlusOne;

template<>
InputParameters validParams<ValuePlusOne>();

class ValuePlusOne : public ChangeVariableData
{
public:
  ValuePlusOne(const std::string & name, InputParameters parameters);
  virtual void initialize();
  virtual void modifySolutionVector();
  virtual Real getValue();

protected:
private:
  Real _foo;

};

#endif //VALUEPLUSONE_H
