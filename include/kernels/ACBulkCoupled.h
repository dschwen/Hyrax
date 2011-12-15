/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  16 November 2011
*
*  This code inherits from CHBulk in ELK
*  
*  This code handles the materials parameters for a coupled 
*  conserved order parameter, non-conserved order parameter
*  system. 
*************************************************************************/

#ifndef ACBulkCoupled_H
#define ACBulkCoupled_H

#include "ACBulk.h"

/* This kernel couples the bulk Alan-Cahn equation term with order parameter eta to the conserved field 
variable term from the Cahn-Hilliard equation (typically concentration)*/


//Forward Declarations
class ACBulkCoupled;

template<>
InputParameters validParams<ACBulkCoupled>();

class ACBulkCoupled : public ACBulk
{
public:

  ACBulkCoupled(const std::string & name, InputParameters parameters);

protected:

  std::string _a2_name;
  std::string _a3_name;
  std::string _a4_name;
  std::string _c2_name;

  virtual Real computeDFDOP(PFFunctionType type);

private:

  MaterialProperty<Real> & _a2;
  MaterialProperty<Real> & _a3;
  MaterialProperty<Real> & _a4;
  MaterialProperty<Real> & _c2;

  VariableValue & _coupled_var;
};

#endif //ACBulkCoupled_H

