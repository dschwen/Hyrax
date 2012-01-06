/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  22 December 2011
*
*************************************************************************/

#ifndef DIRACFORCEDAMR_H
#define DIRACFORCEDAMR_H

#include "ConstantPointSource.h"

#include "mesh_refinement.h"

class DiracForcedAMR;

template<>
InputParameters validParams<DiracForcedAMR>();

/**
 * DiracForcedAMR is designed to have the introduction of a point source after some amount of time
 * has passed.  This kernel is designed to work with the TransientMultiAMR executioner for multiple
 * refinement levels in one timestep.
 */

class DiracForcedAMR : public ConstantPointSource
{
public:
  DiracForcedAMR(const std::string & name, InputParameters parameters);

  virtual void addPoints();
  virtual Real computeQpResidual();

protected:

private:
  int _active_after;
  int _active_for;
};

#endif //DIRACFORCEDAMR_H
