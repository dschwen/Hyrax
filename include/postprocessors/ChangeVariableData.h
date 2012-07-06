#ifndef CHANGEVARIABLEDATA_H
#define CHANGEVARIABLEDATA_H

#include "GeneralPostprocessor.h"

//Forward Declarations
class ChangeVariableData;
class MooseVariable;

template<>
InputParameters validParams<ChangeVariableData>();

class ChangeVariableData : public GeneralPostprocessor
{
public:
  ChangeVariableData(const std::string & name, InputParameters parameters);

  virtual void initialize();
  virtual void execute();

  virtual Real getValue();
  virtual void threadJoin(const Postprocessor & y);

private:
  /// The reference to the variable we want to retrieve in the solution vector
  MooseVariable & _moose_variable;

  /// The number of the variable in the Nonlinear System
  unsigned int _variable_number;
};

#endif // CHANGEVARIABLEDATA_H
