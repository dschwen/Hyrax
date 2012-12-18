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

template<>
InputParameters validParams<MeshSolutionModify>();

class MeshSolutionModify: public Transient
{
public:

  MeshSolutionModify(const std::string & name, InputParameters parameters);

  virtual void endStep();

protected:
  unsigned int _adapt_cycles;
  unsigned int _max_h_level;
};

#endif //MESHSOLUTIONMODIFY_H
