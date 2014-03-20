# This input file demonstrates the coupled multiple Allen-Cahn, Cahn-Hilliard
# equations and explicit nucleation.  It tests calculation of nucleation rate
# with mesh adaptivity.

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 2
  ny = 2
  nz = 0
  xmin = 0
  xmax = 8
  ymin = 0
  ymax = 8
  zmin = 0
  zmax = 0
  elem_type = QUAD4
 uniform_refine = 1
[]

[Variables]
  [./concentration]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
  [../]

  [./n1]
    order = FIRST
    family = LAGRANGE
  [../]
 []

[AuxVariables]
  [./nodal_Supersaturation]
    order = FIRST
    family = LAGRANGE
  [../]

  [./nodal_NucleationRate]
    order = FIRST
    family = LAGRANGE
  [../]

[Kernels]
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
    variable = nodal_Supersaturation
    coupled_var = concentration
    functional_c1 = 0.006
    execute_on = timestep
  [../]

  [./NucleationRate]
    type = AuxNucleationRate
    variable = nodal_NucleationRate
    coupled_aux_var = nodal_Supersaturation
    #Kn1 = 0.008
    #Kn2 = 0.3

    gamma = 0.18
    scale_factor = 900e-22

    Z = 0.1
    Beta_star = 100
    linear_density = 5
    n_OP_vars = 1
    OP_var_names = 'n1'
    execute_on = timestep
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
  type = Steady

  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'


  print_linear_residuals = true

[]

[Outputs]
  file_base = AuxNucleationRate_refine1
  output_initial = true
  exodus = true
  [./console]
    type = Console
    perf_log = true
    linear_residuals = true
  [../]
[]
