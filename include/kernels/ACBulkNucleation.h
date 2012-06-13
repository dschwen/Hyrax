/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  12 June 2012
*
*************************************************************************/

#ifndef ACBULKNUCLEATION_H
#define ACBULKNUCLEATION_H

#include "ACBulkCoupled.h"

//forward declarations
class ACBulkNucleation;

template<>
InputParameters validParams<ACBulkNucleation>();

class ACBulkNucleation : public ACBulkCoupled
{
public:

  ACBulkNucleation(const std::string & name, InputParameters parameters);

protected:

  virtual Real computeDFDOP(PFFunctionType type);

private:

  Real _start_time;
  Real _end_time;
  Real _radius;
//  Real _int_width;
  Point _nucleation_center;
};

#endif //ACBULKNUCLEATION_H
