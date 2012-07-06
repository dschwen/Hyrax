#ifndef CHANGEVARIABLEDATA_H
#define CHANGEVARIABLEDATA_H

#include "GeneralPostprocessor.h"

//Forward Declarations
class ChangeVariableData;
class MooseMesh;
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
  /// A reference to the mesh
  MooseMesh & _mesh;
  /// The reference to the variable we want to retrieve in the solution vector
  MooseVariable & _moose_variable;
  /// The number of the variable in the Nonlinear System
  unsigned int _variable_number;

  unsigned int _foo;
};

#endif // CHANGEVARIABLEDATA_H
