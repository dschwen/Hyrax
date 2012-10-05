[Mesh]
# specify if generated mesh, mesh generated from Cubit, etc.
  type = GeneratedMesh
  dim = 2
  nx = 30
  ny = 30
  nz = 0
  xmin = 0
  xmax = 38.4
  ymin = 0
  ymax = 38.4
  zmin = 0
  zmax = 0
  elem_type = QUAD4
[]

[Variables]
  [./n1]
    order = FIRST
    family = LAGRANGE
  [../]

  [./conc]
    order = THIRD
    family = HERMITE
  [../]
[]

[AuxVariables]
  [./n1_untouched]
    order = FIRST
    family = LAGRANGE
  [../]
[../]

[Kernels]
  [./transient]
   type = TimeDerivative
   variable = n1
  [../]

  [./transient2]
   type = TimeDerivative
   variable = conc
  [../]

  [./CHSolid]
    type = CHBulkCoupled
    variable = conc
    mob_name = M
    coupled_OP_var = n1
  [../]

  [./CHInterface]
    type = CHInterface
    variable = conc
    kappa_name = kappa_c
    mob_name = M
    grad_mob_name = grad_M
  [../]

  [./ACSolid]
    type = ACBulkCoupled
    variable = n1
    coupled_CH_var = conc
    mob_name = L
  [../]

  [./ACInterface]
    type = ACInterface
    variable = n1
    mob_name = L
    kappa_name = kappa_n
  [../]
[]

[AuxKernels]
  [./Report_n1]
   type = ReporterAux
   variable = n1_untouched
   coupled = n1
  execute_on = timestep
  []
[../]

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

[Postprocessors]
  [./Nucleation]
    type = OneSeed
    variables = 'n1'
    radius = 1.8
    seed_value = 1.6
    int_width = 0.9
   # x_position = 38.4
   # y_position = 38.4
    x_position = 19.2
    y_position = 19.2
    restart = true
    restart_time = 0.03
    #output = file
    dwell_time = 0.0011
  [../]
[]

[Executioner]
# how MOOSE should be executed is specified here.
  type = Transient
  scheme = 'bdf2'
#  petsc_options = '-snes_mf_operator -ksp_monitor'
  petsc_options = '-snes_mf_operator'
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = 'hypre boomeramg 101'

  l_max_its = 30
  nl_max_its = 15
  #nl_abs_tol = 1.1e-5

  start_time = 0.0
  num_steps = 5
  dt = 1.0e-2
  abort_on_solve_fail = true

#  [./Adaptivity]
#   coarsen_fraction = 0.05
#   refine_fraction = 0.1
#   max_h_level = 3
#  []
[]

[Output]
# what output should come out of the application is specified here.
  file_base = OneSeed_restart_out
  output_initial = true
  interval = 1
  exodus = true
  perf_log = true
  #postprocessor_csv = true
[]
