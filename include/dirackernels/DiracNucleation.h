/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  15 December 2011
*
*************************************************************************/

#ifndef DIRACNUCLEATION_H
#define DIRACNUCLEATION_H

// Moose Includes
#include "DiracKernel.h"

// Libmesh Includes
#include "mesh_refinement.h"

/* Remember, in the input file, tell this Dirac Kernel that it's operating on the appropriate
   variable - in this case, on the order parameter variable.  When I've got multiple OPs, this is
   going to get interesting. */

//Forward Declarations
class DiracNucleation;

template<>
InputParameters validParams<DiracNucleation>();

/**
 *  DiracNucleation works with the AuxNucleation etc. system to determine where and when Dirac delta "spikes"
 *  are introduced into the order parameter field variable.  This is for the simulation of the explicit
 *  introduction of nuclei as delta functions for the concurrent nucleation and growth algorithm first 
 *  proposed by J.P. Simmons (2000). 
 */

class DiracNucleation : public DiracKernel
{
public:
  DiracNucleation(const std::string & name, InputParameters parameters);

  /**
   * addPoints()
   * @return void function that adds Dirac spikes where  the nucleation field variable says nucleation
   * occured during that time step.
   */
  virtual void addPoints();
  virtual Real computeQpResidual();

protected:
  Real _value;  ///< input value for what you want the Dirac spike to be (probably something like 10.0)
// Maybe some other things here

private: 
  VariableValue & _coupled_nucleation;  ///< AuxVariable for the nucleation yes/no field.

//  Elem *_true_element;
//  Point _true_point;
  MeshRefinement _my_refine;

};


#endif //DIRACNUCLEATION_H
