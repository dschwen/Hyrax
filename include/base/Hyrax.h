/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  16 November 2011
*
*
*  This code handles the materials parameters for a coupled
*  conserved order parameter, non-conserved order parameter
*  system.
*************************************************************************/

#ifndef HYRAX_H
#define HYRAX_H

class Factory;
class ActionFactory;
class Syntax;

namespace Hyrax
{
  /**
   * Register this application and any it depends on.
   */
  void registerApps();
  /**
   * Registers all Kernels and BCs
   */
  void registerObjects(Factory & factory);
}

#endif //HYRAX_H
