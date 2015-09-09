[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 10
  ny = 10
  nz = 0
  xmin = 0
  xmax = 10
  ymin = 0
  ymax = 10
  zmin = 0
  zmax = 0
  elem_type = QUAD4

#  uniform_refine = 1
[]

[Variables]
  active = 'concentration n1'

  [./concentration]
    order = THIRD
    family = HERMITE
#    [./InitialCondition]
#      int_width = 1
#      invalue = 0.01
#      outvalue = 0.01
#      radius = 1
#      type = SmoothCircleIC
#      var_name = concentration
#      x1 = 5
#      y1 = 5
#    [../]
  [../]

  [./n1]
    order = FIRST
    family = LAGRANGE
#    [./InitialCondition]
#      int_width = 1
#      invalue = 1.0
#      outvalue = 0.0
#      radius = 1
#      type = SmoothCircleIC
#      var_name = n1
#      x1 = 5
#      y1 = 5
#    [../]
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

[DiracKernels]
  active = 'dirac'
  [./dirac]
    type = ConstantPointSource
    variable = n1
    value = 10.0
    point = '6 6 0'
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
  type = TransientMultiAMR
  scheme = 'crank-nicolson'

  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'




  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = 'hypre boomeramg 101'

  l_max_its = 15
  nl_max_its = 10

  start_time = 0.0
  num_steps = 5
  dt = 0.03

  num_refines = 2

  [./Adaptivity]
    error_estimator = LaplacianErrorEstimator
    refine_fraction = 0.98
    coarsen_fraction = 0.0
    max_h_level = 6
  [../]
[]

[Outputs]
  file_base = out
  exodus = true
[]ls
