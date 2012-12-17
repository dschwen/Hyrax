/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  13 December 2012
*
*************************************************************************/

#include "Nucleus.h"

#include "libmesh.h"
#include "point.h"

Nucleus::Nucleus()
{
  _location.zero();
  _start_time = 0.0;
  _end_time = 0.0;
  _orientation = 0;
}

Nucleus::Nucleus(const Nucleus &a)
{
  *this = a;
}

Nucleus &
Nucleus::operator= (const Nucleus &a)
{
  _location = a._location;
  _start_time = a._start_time;
  _end_time = a._end_time;
  _orientation = a._orientation;

  return *this;
}

Point
Nucleus::getLocation() const
{
  return _location;
}

Real
Nucleus::getStartTime() const
{
  return _start_time;
}

Real
Nucleus::getEndTime() const
{
  return _end_time;
}

int
Nucleus::getOrientation() const
{
  return _orientation;
}

void
Nucleus::setLocation(Point a)
{
  _location = a;
}

void
Nucleus::setStartTime(Real a)
{
  _start_time = a;
}

void
Nucleus::setEndTime(Real a)
{
  _end_time = a;
}

void
Nucleus::setOrientation(int a)
{
  _orientation = a;
}



