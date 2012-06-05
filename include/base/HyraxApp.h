/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  30 May 2012
*
*  
*  This code handles the materials parameters for a coupled 
*  conserved order parameter, non-conserved order parameter
*  system. 
*************************************************************************/

#ifndef HYRAXAPP_H
#define HYRAXAPP_H

#include "MooseApp.h"

class HyraxApp : public MooseApp
{
public:
  HyraxApp(int argc, char * argv []);
};

#endif //HYRAXAPP_H
