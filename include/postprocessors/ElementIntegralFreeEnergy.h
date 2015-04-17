/****************************************************************/
/****************************************************************/

#ifndef ELEMENTINTEGRALFREEENERGY_H
#define ELEMENTINTEGRALFREEENERGY_H

#include "ElementIntegralPostprocessor.h"
//#include "MooseVariableInterface.h"

class ElementIntegralFreeEnergy;

template<>
InputParameters validParams<ElementIntegralFreeEnergy>();

class ElementIntegralFreeEnergy :
  public ElementIntegralPostprocessor//,
//  public MooseVariableInterface
{
public:
  ElementIntegralFreeEnergy(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeQpIntegral();

  MaterialProperty<Real> & _phase1_energy;
  MaterialProperty<Real> & _phase2_energy;

  MaterialProperty<Real> & _grad_OP_term;
  MaterialProperty<Real> & _grad_C_term;

  MaterialProperty<Real> & _Heaviside;

};

#endif //ELEMENTINTEGRALFREEENERGY_H
