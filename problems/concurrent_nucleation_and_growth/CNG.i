[Mesh]
  [./Generation]
    dim = 2
    nx = 20
    ny = 20
    nz = 0
    xmin = 0
    xmax = 40
    ymin = 0
    ymax = 40
    zmin = 0
    zmax = 0
    elem_type = QUAD4
  [../]
  uniform_refine = 2
[]

[Variables]
  active = 'concentration n1'

  [./concentration]
    order = THIRD
    family = HERMITE
    [./InitialCondition]
      value = 0.1
      type = ConstantIC
      var_name = concentration
    [../]
  [../]

  [./n1]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = ConstantIC
      var_name = n1
      value = 0.0
    [../]
  [../]

[]

[AuxVariables]
  active = 'nucleation supersaturation j_star p_nm'

  [./nucleation]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./supersaturation]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./j_star]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./p_nm]
    order = CONSTANT
    family = MONOMIAL
  [../]

[]

[Kernels]
  active = 'dcdt dn1dt CHSolid CHInterface ACSolid ACInterface'

  [./dcdt]
    type = TimeDerivative
    variable = concentration
  [../]

  [./dn1dt]
    type = TimeDerivative
    variable = n1
  [../]

  [./CHSolid]
    type = CHBulkCoupled
    variable = concentration
    mob_name = M
    coupled_var = n1
    first_landau = A1
    second_landau = A2
    first_well = C1
  [../]

  [./CHInterface]
    type = CHInterface
    variable = concentration
    kappa_name = kappa_c
    mob_name = M
    grad_mob_name = grad_M
  [../]

  [./ACSolid]
    type = ACBulkCoupled
    variable = n1
    mob_name = L
    coupled_var = concentration
    second_landau = A2
    third_landau = A3
    fourth_landau = A4
    second_well = C2
  [../]

  [./ACInterface]
    type = ACInterface
    variable = n1
    mob_name = L
    kappa_name = kappa_n
  [../]
[]

[AuxKernels]
  active = 'calc_nuc calc_supersat calc_j_star calc_p_nm'

  [./calc_nuc]
    type = AuxNucleation
    variable = nucleation
    coupled_aux_var = p_nm
  [../]

  [./calc_supersat]
    type = AuxSupersaturation
    variable = supersaturation
    coupled_var = concentration
  [../]

  [./calc_j_star]
    type = AuxNucleationRate
    variable = j_star
    Kn1 = 0.001
    Kn2 = 0.33
    coupled_aux_var = supersaturation
  [../]

  [./calc_p_nm]
    type = AuxNucleationProbability
    variable = p_nm
    coupled_aux_var = j_star
  [../]
[]

[DiracKernels]
  active = 'dirac'
  [./dirac]
    type = DiracNucleation
    variable = n1
    value = 5.0
    nucleation = nucleation
  [../]
[]

[BCs]
active = 'Periodic'
  [./Periodic]
    [./left_right]
      primary = 0
      secondary = 2
      translation = '0 10 0'
    [../]

    [./top_bottom]
      primary = 1
      secondary = 3
      translation = '-10 0 0'
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
  petsc_options = '-snes_mf_operator -ksp_monitor'

  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = 'hypre boomeramg 101'

  l_max_its = 15
  nl_max_its = 10

#  e_max = 6.0e2
#  e_tol = 5.0e-1

  start_time = 0.0
  num_steps = 15
  dt = 0.0003

# The adaptivity section should probably be played with at will.
#  [./Adaptivity]
#    error_estimator = LaplacianErrorEstimator
#    refine_fraction = 0.85
#    coarsen_fraction = 0.1
#    max_h_level = 4
#  [../]
[]

[Output]
  file_base = out
  output_initial = true
  interval = 1
  exodus = true
  perf_log = true
[]


