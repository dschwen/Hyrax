# this input file is to test the concurrent nucleation and growth with mesh adaptivity

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 10
  ny = 10
  nz = 0
  xmin = 0
  xmax = 76.8 #0.3*256
  ymin = 0
  ymax = 76.8
  zmin = 0
  zmax = 0
  elem_type = QUAD4
  #uniform_refine = 3 # 80 elements, dx=0.96...
[]

[Variables]
  [./concentration]
    order = THIRD
    family = HERMITE
      [./InitialCondition]
      type = ConstantIC
      value = 0.1
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
    #order = FIRST
    #family = LAGRANGE
  [../]

 [./elemental_NucleationProbability]
    order = CONSTANT
    family = MONOMIAL
    #order = FIRST
    #family = LAGRANGE
  [../]

 [./elemental_NucleationRate]
    order = CONSTANT
    family = MONOMIAL
    #order = FIRST
    #family = LAGRANGE
  [../]
[]

[Kernels]
  [./dcdt]
    type = TimeDerivative
    variable = concentration
  [../]

  [./dn1dt]
    type = TimeDerivative
    variable = n1
  [../]

  [./CHSolid]
   # type = CHBulkPolyCoupled
    type = CHBulkCoupled
    variable = concentration
    mob_name = M
   # n_OP_variables = 1
   # OP_variable_names = 'n1'
    coupled_OP_var = n1
  [../]

  [./CHInterface]
    type = CHInterface
    variable = concentration
    kappa_name = kappa_c
    mob_name = M
    grad_mob_name = grad_M
  [../]

  [./ACSolidn1]
  #  type = ACBulkPolyCoupled
    type = ACBulkCoupled
    variable = n1
    mob_name = L
    coupled_CH_var = concentration
   # n_OP_vars = 1
   # OP_var_names = 'n1
   # OP_number = 1
  [../]

  [./ACInterfacen1]
    type = ACInterface
    variable = n1
    mob_name = L
    kappa_name = kappa_n
  [../]
[]

[AuxKernels]
  [./Supersaturation]
    type = AuxSupersaturation
    variable = elemental_Supersaturation
    coupled_var = concentration
    functional_c1 = 0.006
    execute_on = timestep
  [../]

  [./NucleationRate]
    type = AuxNucleationRate
    variable = elemental_NucleationRate
    OP_var_names = 'n1'
    n_OP_vars = 1
    coupled_aux_var = elemental_Supersaturation
    Beta_star = 0.1
    linear_density = 70
    Z = 0.1
    Kn2 = 0.3
    execute_on = timestep
  [../]

  [./NucleationProbability]
    type = AuxNucleationProbability
    variable = elemental_NucleationProbability
    coupled_aux_var = elemental_NucleationRate
    coupled_variables = 'n1'
    n_OP_vars = 1
    execute_on = timestep
  [../]
[]

[BCs]
  [./c_BC]
    type = NeumannBC
    variable = concentration
    boundary = '0 1 2 3'
    value = 0.0
  [../]

  [./n1_BC]
    type = NeumannBC
    variable = n1
    boundary = '0 1 2 3'
    value = 0.0
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

[UserObjects]
  [./NLUO]
    type = NucleationLocationUserObject
    variable = n1
    coupled_aux = elemental_NucleationProbability
    dwell_time = 0.1
    num_orientations = 1
  [../]

  [./NISM]
    type = NucleusIntroductionSolutionModifier
    nucleation_userobject = NLUO
    variables = 'n1'
    num_vars = 1
    seed_value = 1.6
    radius = 1.8
    int_width = 0.9
    execute_on = custom
  [../]
[../]

[Executioner]
  type = MeshSolutionModify
  scheme = 'crank-nicolson'
  num_steps = 10
  dt = 0.01
  adapt_cycles = 2
#  max_h_level = 2

  petsc_options = -snes_mf_operator
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'

#  [./Adaptivity]
#   coarsen_fraction = 0.01
#   refine_fraction = 0.8
#   max_h_level = 5
#   error_estimator = LaplacianErrorEstimator
#  [../]
[../]

[Adaptivity]
   marker = NM
 [./Markers]
    [./NM]
      type = NucleationMarker
      nucleation_userobject = NLUO
    [../]
    [./EFM]
      type = ErrorFractionMarker
      coarsen = 0.01
      refine = 0.8
      indicator = GJI
    [../]
    [./combo]
      type = ComboMarker
      markers = 'NM EFM'
    [../]
 [../]

 [./Indicators]
   [./GJI]
    type = GradientJumpIndicator
    variable = concentration
   [../]
 [../]
[]

[Output]
  file_base = testadapt
  output_initial = true
  interval = 1
  exodus = true
  perf_log = true
[../]