/****************************************************************/
/* Stress BC:                                                   */
/*     This BC is intended only for testing purpose.            */
/*     Mathematically, Stress BC is not well-posed.             */
/*     Anyway, it reads six functions:                          */
/*     stress_xx, stress_xy, stress_yy,                         */
/*     stress_yz, stress_zx, stress_zz                          */
/*     Then, it computes the traction by multiplying            */
/*     stress tensor with the normal to get the traction.       */
/*     The traction is used for the problem.                    */
/*                                                              */
/****************************************************************/

#ifndef STRESSBC_H
#define STRESSBC_H

#include "IntegratedBC.h"

//Forward Declarations
class StressBC;

template<>
InputParameters validParams<StressBC>();

class StressBC : public IntegratedBC
{
public:

  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  StressBC(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeQpResidual();

private:
  int _component; 
  Function & _stress_xx;
  Function & _stress_xy;
  Function & _stress_yy;
  Function & _stress_yz;
  Function & _stress_zx;
  Function & _stress_zz;
};

#endif
