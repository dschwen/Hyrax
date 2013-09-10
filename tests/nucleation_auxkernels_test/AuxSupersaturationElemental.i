# This input file demonstrates the AuxSupersaturation aux kernel with elemental auxvars.

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
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      variable = concentration
    [../]
  [../]
[]

[AuxVariables]
  [./elemental_Supersaturation]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Kernels]
  [./Diffusion]
    type = Diffusion
    variable = concentration
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
[]

[BCs]
  [./c_BC]
    type = DirichletBC
    variable = concentration
    boundary = '0 1'
    value = 1.0
  [../]

  [./c_BC2]
    type = DirichletBC
    variable = concentration
    boundary = '2 3'
    value = 0.0
  [../]
[]


[Executioner]
  type = Steady
  #scheme = 'crank-nicolson'

  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'


  print_linear_residuals = true


  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = 'hypre boomeramg 101'

  l_max_its = 15
  nl_max_its = 10
[]

[Output]
  file_base = AuxSupersaturationElemental
  output_initial = true
  exodus = true
  perf_log = true
[]

