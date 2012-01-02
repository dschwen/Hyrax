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
 * DiracForcedAMR is designed to have the introduction of a point source and force mesh adaptivity
 * at the point source. This kernel is designed to work with the TransientMultiAMR executioner. 
 * The mesh refinement for the element in which the point at which the Dirac kernel was introduced 
 * needs to be forced several levels.
 */

class DiracForcedAMR : public ConstantPointSource
{
public:
  DiracForcedAMR(const std::string & name, InputParameters parameters);

  virtual void addPoints();
  virtual Real computeQpResidual();

protected:

private:
//  MeshRefinement _my_refine;
  Real _active_after;
};

#endif //DIRACFORCEDAMR_H
