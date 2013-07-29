/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  24 July 2013
*
*************************************************************************/

#include "OPVariantKernelAction.h"
#include "Factory.h"

#include <iostream>

template<>
InputParameters validParams<OPVariantKernelAction>()
{
  InputParameters params = validParams<Action>();

  params.addRequiredParam<unsigned int>("number_OPs", "specifies the number of order parameters");
  params.addRequiredParam<std::string>("OP_name_base", "name base of the OP variable");
  params.addParam<std::string>("kappa_name_OP", "kappa_n", "the OP kappa used with ACInterface");
  params.addRequiredParam<std::string>("coupled_CH_var", "coupled conserved variable for ACBulkPolyCoupled");

  return params;
}

OPVariantKernelAction::OPVariantKernelAction(const std::string & name, InputParameters params) :
    Action(name, params),
    _num_OPs(getParam<unsigned int>("number_OPs")),
    _OP_name_base(getParam<std::string>("OP_name_base")),
    _kappa_name_OP(getParam<std::string>("kappa_name_OP")),
    _coupled_CH_var(getParam<std::string>("coupled_CH_var"))
{
}

void
OPVariantKernelAction::act()
{
  for(unsigned int i = 1; i <= _num_OPs; ++i)
  {
    //create variable names
    std::string variable_name = _OP_name_base;
    std::stringstream out;
    out << i;
    variable_name.append(out.str());

    //create vector of coupled variable names for each kernel
    std::vector<std::string> coupled_var_name_vector;
    coupled_var_name_vector.resize(_num_OPs - 1);
    //counter for vector
    unsigned int k = 0;

    for (unsigned int j = 1; j <= _num_OPs; ++j)
    {
      if( j != i)
      {
        std::string coupled_var_name = _OP_name_base;
        std::stringstream out2;
        out2 << j;
        coupled_var_name.append(out2.str());
        coupled_var_name_vector[k] = coupled_var_name;
        k++;
      }
    }

    std::cout<<"kernel action before setting params for ACBulkPolyCoupled "<< i <<std::endl;

    std::string kernel_name;

    //Get the parameters for the ACBulkPolyCoupled kernel and kernels to problem
    InputParameters action_params = _factory.getValidParams("ACBulkPolyCoupled");
    action_params.set<NonlinearVariableName>("variable") = variable_name;
    action_params.set<std::vector<std::string> >("OP_var_names") = coupled_var_name_vector;
    action_params.set<std::string>("coupled_CH_var") = _coupled_CH_var;
    action_params.set<int>("n_OP_vars") = _num_OPs;
    action_params.set<int>("OP_number") = i;

    kernel_name = "ACBulkPolyCoupled_";
    kernel_name.append(variable_name);

    _problem->addKernel("ACBulkPolyCoupled", kernel_name, action_params);

    std::cout<<"kernel action after setting params for ACBulkPolyCoupled"<<std::endl;

    //Get the parameters for the ACInterface kernel and add kernels to problem
    action_params = _factory.getValidParams("ACInterface");
    action_params.set<NonlinearVariableName>("variable") = variable_name;
    action_params.set<std::string>("kappa_name") = _kappa_name_OP;

    kernel_name = "ACInterface_";
    kernel_name.append(variable_name);
    _problem->addKernel("ACInterface", kernel_name, action_params);

    //Get the parameters for the TimeDerivative kernel and add kernels to problem
    action_params = _factory.getValidParams("TimeDerivative");
    action_params.set<NonlinearVariableName>("variable") = variable_name;

    kernel_name = "TimeDeriv_";
    kernel_name.append(variable_name);

    _problem->addKernel("TimeDerivative", kernel_name, action_params);


    //anything else you need to add?
  }
}
