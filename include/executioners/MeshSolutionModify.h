/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  11 December 2012
*
*  This code is originally by Derek Gaston, in moose_test - just been
*  copied and minorly modified.
*
*************************************************************************/

#ifndef MESHSOLUTIONMODIFY_H
#define MESHSOLUTIONMODIFY_H

//#include "Transient.h"
#include "SolutionTimeAdaptive.h"

// Forward Declarations
class MeshSolutionModify;
class NucleationLocationUserObject;

template<>
InputParameters validParams<MeshSolutionModify>();

class MeshSolutionModify: public SolutionTimeAdaptive //Transient
{
public:

  MeshSolutionModify(const std::string & name, InputParameters parameters);

//  virtual void takeStep(Real input_dt = -1.0);
  virtual void endStep();
  virtual void preExecute();
  Real computeDT();

protected:
  unsigned int _adapt_cycles;
  unsigned int _adapt_nucleus;

  bool _use_nucleation_userobject;
  bool _new_nucleus;

  const NucleationLocationUserObject*  _nucleation_userobject;
};

#endif //MESHSOLUTIONMODIFY_H
