# This input file is designed to test the AC*Nucleation kernels along with the ValueNucleation
# kernel and the ForcingFunctionNucleation kernel, for performing an L2 projection to make a
# nucleus.

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
  [./n1]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = ConstantIC
      value = 0.01
    [../]
  [../]

  [./concentration]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = ConstantIC
      value = 0.1
    [../]
  [../]
[]

[Functions]
  [./forcing_function]
    type = ParsedFunction
    value = 1.6
  [../]
[]

[Kernels]
  [./dn1dt]
    type = TimeDerivative
    variable = n1
  [../]

  [./ACSolid]
    type = ACBulkNucleation
    variable = n1
    coupled_CH_var = concentration
    mob_name = L
    center_point = '50 50 0'
    start_time = 0.0011
    end_time = 0.02
    radius = 10.0
  [../]

  [./ACInterface]
    type = ACInterfaceNucleation
    variable = n1
    mob_name = L
    kappa_name = kappa_n
    center_point = '50 50 0'
    start_time = 0.0011
    end_time = 0.02
    radius = 10.0
  [../]

  [./ValueNucleation_test]
    type = ValueNucleation
    variable = n1
    center_point = '50 50 0'
    radius = 10.0
  #  int_width = 5.0
    start_time = 0.0011
    end_time = 0.02
  [../]

  [./ForcingNucleation_test]
    type = ForcingFunctionNucleation
    variable = n1
    function = forcing_function
    center_point = '50 50 0'
    radius = 10.0
 #   int_width = 5.0
    start_time = 0.0011
    end_time = 0.02
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
  file_base = testACNucleationForcing
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
