# This input file demonstrates the AuxSupersaturation aux kernel. It's NODAL - this is going to
# to have to change soon!

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 10
  ny = 10
  nz = 0
  xmin = 0
  xmax = 50
  ymin = 0
  ymax = 50
  zmin = 0
  zmax = 0
  elem_type = QUAD4
[]

[Variables]
  [./concentration]
    order = THIRD
    family = HERMITE
    [./InitialCondition]
      type = RandomIC
      variable = concentration
    [../]
  [../]

 [./n1]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      variable = n1
    [../]
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

  [./nodal_NucleationProbability]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
  [./CH]
    type = CHBulkCoupled
    variable = concentration
    mob_name = M
    coupled_OP_var = n1
  [../]

  [./AC]
    type = ACBulkCoupled
    variable = n1
    coupled_CH_var = concentration
    mob_name = L
  [../]

  [./dcdt]
    type = TimeDerivative
    variable = concentration
  [../]

  [./dndt]
    type = TimeDerivative
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
    OP_var_names = 'n1'
    n_OP_vars = 1
    execute_on = timestep
  [../]

 [./AuxNucleationProbability]
    type = AuxNucleationProbability
    variable = nodal_NucleationProbability
    coupled_aux_var = nodal_NucleationRate
    coupled_variables = 'n1'
    n_OP_vars = 1
    execute_on = timestep
 [../]
[]

[BCs]
  [./Periodic]
    [./all]
      auto_direction = 'x y'
    [../]
  [../]
[]

[Materials]
  [./constant]
    type = PFMobilityLandau
    block = 0
    mob_CH = 0.4
    mob_AC = 0.4
    kappa_CH = 1.5
    kappa_AC = 1.5
    A1 = 18.5
    A2 = -8.5
    A3 = 11.5
    A4 = 4.5
    C1 = 0.006
    C2 = 0.59
  [../]
[]


[Executioner]
  type = Transient
  scheme = 'crank-nicolson'

  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'

  petsc_options = '-ksp_monitor'

  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = 'hypre boomeramg 101'

  l_max_its = 15
  nl_max_its = 10

  start_time = 0.0
  num_steps = 2
  dt = 0.001
[]

[Output]
  file_base = AuxNucleationProbability
  output_initial = true
  exodus = true
  perf_log = true
[]

