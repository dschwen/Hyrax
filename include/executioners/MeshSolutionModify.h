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

#include "Transient.h"

// Forward Declarations
class MeshSolutionModify;
class NucleationLocationUserObject;

template<>
InputParameters validParams<MeshSolutionModify>();

class MeshSolutionModify: public Transient
{
public:

  MeshSolutionModify(const std::string & name, InputParameters parameters);

  virtual void endStep();
  //virtual void preExecute();

protected:
  unsigned int _adapt_cycles;
  unsigned int _adapt_nucleus;

  //const NucleationLocationUserObject*  _nucleation_userobject;
};

#endif //MESHSOLUTIONMODIFY_H
