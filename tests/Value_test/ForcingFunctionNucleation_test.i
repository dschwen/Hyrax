# This input file is designed to test the Value kernel using the
#UserForcingFunction kernel, which will is designed
# to be used for performing an L2 projection.

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 25
  ny = 25
  nz = 0
  xmin = 0
  xmax = 100
  ymin = 0
  ymax = 100
  zmin = 0
  zmax = 0
  elem_type = QUAD4
[]

[Variables]
  [./diffused]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = SmoothCircleIC
      int_width = 1.5
      invalue = 0.6
      outvalue = 0.1
      radius = 5.0
      x1 = 50
      y1 = 50
      z1 = 0
    [../]
  [../]
[]

[Functions]
  [./forcing_function]
    type = ParsedFunction
    value = 1.0*x
  [../]
[]

[Kernels]
  [./Value_test]
    type = Value
    variable = diffused
  [../]

  [./Forcing_test]
    type = UserForcingFunction
    variable = diffused
    function = forcing_function
  [../]
[]

[BCs]
  [./Periodic]
    [./all]
     auto_direction = 'x y'
    [../]
  [../]
[]

[Executioner]
  type = Steady

  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'


[]

[Outputs]
  file_base = testValueForcing
  output_initial = true
  exodus = true
  [./console]
    type = Console
    perf_log = true
    linear_residuals = true
  [../]
[]
