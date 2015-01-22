# This input file demonstrates the coupled multiple Allen-Cahn, Cahn-Hilliard
# equations and explicit nucleation.  It tests calculation of nucleation rate
# with mesh adaptivity.  Uses ELEMENTAL aux variables.

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 4
  ny = 4
  nz = 0
  xmin = 0
  xmax = 8
  ymin = 0
  ymax = 8
  zmin = 0
  zmax = 0
  elem_type = QUAD4
#  uniform_refine = 1
[]

[Variables]
  [./concentration]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
    #  value = 1.0
    [../]
  [../]

  [./n1]
   order = FIRST
   family = LAGRANGE
   [./InitialCondition]
      type = RandomIC
      min = 0.0
      max = 0.2
    [../]
  [../]
 []

[AuxVariables]
  [./elemental_Supersaturation]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./elemental_NucleationRate]
    order = CONSTANT
    family = MONOMIAL
  [../]

[Kernels]
  [./time_deriv_diff]
    type = TimeDerivative
    variable = concentration
  [../]

  [./time_deriv_diff2]
    type = TimeDerivative
    variable = n1
  [../]

 [./diffusion]
    type = Diffusion
    variable = concentration
  [../]

  [./diff2]
   type = Diffusion
   variable = n1
  [../]
[]

[AuxKernels]
  [./Supersaturation]
    type = AuxSupersaturation
    variable = elemental_Supersaturation
    coupled_var = concentration
    functional_c1 = 0.006
    execute_on = timestep_end
  [../]

  [./NucleationRate]
    type = AuxNucleationRate
    variable = elemental_NucleationRate
    coupled_aux_var = elemental_Supersaturation
    #Kn1 = 0.008
    #Kn2 = 0.3

    gamma = 0.18
    scale_factor = 900e-22

    Z = 0.1
    Beta_star = 100
    linear_density = 5
    n_OP_vars = 1
    OP_var_names = 'n1'
    execute_on = timestep_end
  [../]
[]

[BCs]
 [./Periodic]
  [./all]
    variable = concentration
    auto_direction = 'x y'
  [../]
 [../]
[]

[Executioner]
  type = Transient
   dt = 0.001
   num_steps = 2

  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'


  print_linear_residuals = true

[]

[Outputs]
  file_base = AuxNucleationRateElemental_refine0
  output_initial = true
  exodus = true
  print_linear_residuals = true
  print_perf_log = true
[]
