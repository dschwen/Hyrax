#ifndef CHANGEVARIABLEDATA_H
#define CHANGEVARIABLEDATA_H

#include "GeneralPostprocessor.h"

//Forward Declarations
class ChangeVariableData;
class MooseMesh;
class MooseVariable;
class NonlinearSystem;

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
  /// A reference to the nonlinear system
  NonlinearSystem & _nl;
  /// The number of the variable in the Nonlinear System
  unsigned int _variable_number;
  /// This is the storage location for the value we will return in this postprocessor
  Real _foo;
  /// A reference to the coupled variable
  MooseVariable & _coupled;
};

#endif // CHANGEVARIABLEDATA_H
