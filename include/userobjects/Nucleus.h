/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  13 December 2012
*
*   Nucleus is designed to hold the information needed for a nucleus in
*   a microstructure.
*
*************************************************************************/

#ifndef NUCLEUS_H
#define NUCLEUS_H

#include "libmesh.h"
#include "point.h"

class Nucleus
{
public:

  /**
   * Default constructor, everything initialized to zero
   */
  Nucleus();

  /**
   * Copy constructor
   */
  Nucleus(const Nucleus &a);

  Nucleus & operator= (const Nucleus &a);

  ~Nucleus() {}

  Point getLocation() const;

  Real getStartTime() const;

  Real getEndTime() const;

  int getOrientation() const;

  void setLocation(Point a);

  void setStartTime(Real a);

  void setEndTime(Real a);

  void setOrientation(int a);

protected:

private:

  Point _location; //holds the central location of the nucleus
  Real _start_time; //when this nucleus came into existence
  Real _end_time; // when the nucleation event ends - NOT necessarily the lifetime of the nucleus
  int _orientation; // the orientation of the particle
};

#endif //NUCLEUS_H
