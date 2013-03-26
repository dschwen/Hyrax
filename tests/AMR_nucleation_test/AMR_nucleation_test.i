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
[]

[Variables]
  [./concentration]
    order = FIRST
    family = LAGRANGE
      [./InitialCondition]
      type = RandomIC
      min = 0.007
      max = 0.01
    [../]
  [../]

#  [./n1]
#    order = FIRST
#    family = LAGRANGE
#    [./InitialCondition]
#      type = RandomIC
#      min = 0.0
#      max = 0.2
#    [../]
#  [../]
[]

[AuxVariables]
  [./elemental_Supersaturation]
    order = CONSTANT
    family = MONOMIAL
   [../]

 [./elemental_NucleationProbability]
    order = CONSTANT
    family = MONOMIAL
  [../]

 [./elemental_NucleationRate]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Kernels]
  [./dcdt]
    type = TimeDerivative
    variable = concentration
  [../]

  [./diff]
    type = Diffusion
    variable = concentration
  [../]
#  [./dn1dt]
#    type = TimeDerivative
#    variable = n1
#  [../]

#  [./CHSolid]
#    type = CHBulkCoupled
#    variable = concentration
#    mob_name = M
#    coupled_OP_var = n1
#  [../]

#  [./CHInterface]
#    type = CHInterface
#    variable = concentration
#    kappa_name = kappa_c
#    mob_name = M
#    grad_mob_name = grad_M
#  [../]

#  [./ACSolidn1]
#    type = ACBulkCoupled
#    variable = n1
#    mob_name = L
#    coupled_CH_var = concentration
#  [../]

#  [./ACInterfacen1]
#    type = ACInterface
#    variable = n1
#    mob_name = L
#    kappa_name = kappa_n
#  [../]
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
    OP_var_names = 'concentration'
    n_OP_vars = 1
    coupled_aux_var = elemental_Supersaturation
    Beta_star = 0.01
    linear_density = 80
    Z = 0.01
    #Kn2 = 0.3

    gamma = 0.18
    Kb = 1
    scale_factor = 1
    execute_on = timestep
  [../]

  [./NucleationProbability]
    type = AuxNucleationProbability
    variable = elemental_NucleationProbability
    coupled_aux_var = elemental_NucleationRate
    coupled_variables = 'concentration'
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

#  [./n1_BC]
#    type = NeumannBC
#    variable = n1
#    boundary = '0 1 2 3'
#    value = 0.0
#  [../]
[]

#[Materials]
#  [./constant]
#    type = PFMobilityLandau
#    block = 0
#    mob_CH = 0.4
#    mob_AC = 0.4
#    kappa_CH = 1.5
#    kappa_AC = 1.5
#    A1 = 18.5
#    A2 = -8.5
#    A3 = 11.5
#    A4 = 4.5
#    C1 = 0.006
#    C2 = 0.59
#  [../]
#[]

[UserObjects]
  [./NLUO]
    type = NucleationLocationUserObject
    #variable = concentration
    coupled_aux_vars = 'elemental_NucleationProbability'
    n_coupled_aux = 1
    dwell_time = 0.1
    num_orientations = 1
  [../]

  [./NISM]
    type = NucleusIntroductionSolutionModifier
    nucleation_userobject = NLUO
    variables = 'concentration'
    num_vars = 1
    seed_value = 1.6
    radius = 1.8
    int_width = 0.9
    execute_on = custom
  [../]
[]

[Executioner]
  type = MeshSolutionModify
  scheme = 'crank-nicolson'
  num_steps = 2
  dt = 0.01
  dtmin = 0.0001
  dtmax = 0.1
  adapt_nucleus = 2

  petsc_options = -snes_mf_operator
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
  use_nucleation_userobject = true
  nucleation_userobject = NLUO
[]

[Adaptivity]
   marker = NM
 [./Markers]
    [./NM]
      type = NucleationMarker
      nucleation_userobject = NLUO
      max_h_level = 2
    [../]
 [../]
[]

[Output]
  file_base = AMR_nucleation_test
  output_initial = true
  interval = 1
  exodus = true
  perf_log = true
[../]
