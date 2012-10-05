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
    [./InitialCondition]
      type = RandomIC
    [../]
  [../]
[]


[Kernels]
  [./transient]
   type = TimeDerivative
   variable = n1
  [../]

  [./diffused]
    type = Diffusion
    variable = n1
  [../]
[]


[BCs]
  [./n1_BC]
    type = NeumannBC
    variable = n1
    boundary = '0 1 2 3'
    value = 0.0
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
    #restart = true
    #restart_time = 0.1
    #output = file
    dwell_time = 0.021
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
  num_steps = 4
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
  file_base = OneSeed_out
  output_initial = true
  interval = 1
  exodus = true
  perf_log = true
  #postprocessor_csv = true
[]
