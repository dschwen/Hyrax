/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  19 December 2012
*
*  ErrorFractionMaxHMarker is a clone of ErrorFractionMarker with the
*  added bonus of setting a maximum refinement level.
*
*************************************************************************/

#ifndef ERRORFRACTIONMAXHMARKER_H
#define ERRORFRACTIONMAXHMARKER_H

#include "ErrorFractionMarker.h"

class ErrorFractionMaxHMarker;

template<>
InputParameters validParams<ErrorFractionMaxHMarker>();

class ErrorFractionMaxHMarker : public ErrorFractionMarker
{
public:
  ErrorFractionMaxHMarker(const std::string & name, InputParameters parameters);
  virtual ~ErrorFractionMaxHMarker(){};

  // virtual void markerSetup();

protected:
  virtual MarkerValue computeElementMarker();

  // Real _coarsen;
  //Real _refine;

  //Real _max;
  //Real _min;
  //Real _delta;
  //Real _refine_cutoff;
  //Real _coarsen_cutoff;
  unsigned int _max_h_level;
};

#endif // ERRORFRACTIONMAXHMARKER_H
