[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 2
  ny = 2
  nz = 0
  zmin = 0
  zmax = 0
  elem_type = QUAD4
[]

[Variables]
  [./u]
    order = FIRST
    family = LAGRANGE
  [../]

  [./v]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
  [./diff]
    type = Diffusion
    variable = u
  [../]

  [./diff2]
    type = Diffusion
    variable = v
  [../]
[]

[BCs]
  [./left_u]
    type = DirichletBC
    variable = u
    boundary = 'left bottom'
    value = 0
  [../]

  [./right_u]
    type = DirichletBC
    variable = u
    boundary = 'right top'
    value = 1
  [../]

  [./left_v]
    type = DirichletBC
    variable = v
    boundary = 'left'
    value = 1
  [../]

  [./right_r]
    type = DirichletBC
    variable = v
    boundary = 'right'
    value = 0
  [../]
[]

[Postprocessors]
  [./change_variable]
    type = ChangeVariableData
    variable = u
    coupled = v
  [../]
[]

[Executioner]
  type = Steady
  petsc_options = '-snes_mf_operator'
[]

[Output]
  output_initial = true
  interval = 1
  exodus = true
[]
