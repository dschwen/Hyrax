# This input file is designed to test the ValueNucleation kernel using the
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
#  uniform_refine = 2
[]

[Variables]
  [./diffused]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = ConstantIC
      value = 0.01
    [../]
  [../]
[]

[Functions]
  [./forcing_function]
    type = ParsedFunction
    value = 100*t
  [../]
[]

[Kernels]
  [./ValueNucleation_test]
    type = ValueNucleation
    variable = diffused
    center_point = '50 50 0'
    radius = 10.0
  #  int_width = 5.0
    start_time = 0.0011
    end_time = 0.02
  [../]

  [./ForcingNucleation_test]
    type = ForcingFunctionNucleation
    variable = diffused
    function = forcing_function
    center_point = '50 50 0'
    radius = 10.0
 #   int_width = 5.0
    start_time = 0.0011
    end_time = 0.02
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
  type = Transient
  petsc_options = '-snes_mf_operator'
  dt = 0.001
  end_time = 0.003
[]


[Output]
  file_base = testValueNucleationForcing
  output_initial = true
  interval = 1
  exodus = true
  perf_log = true
#
#  [./OverSampling]
#   exodus = true
#   refinements = 4
#   output_initial = true
#   interval = 1
#   [../]
[]
