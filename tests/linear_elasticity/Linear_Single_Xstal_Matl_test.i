# This input file is designed to test the LinearSingleCrystalPrecipitate class.  This test is
# for regression testing.  This just takes the material properties and puts them into
# aux variables; the diffusion kernel is just to have a simple kernel to run the test.

# THIS TEST TESTS THE "ELASTICITY_TENSOR" MATERIAL PROPERTY.

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
  uniform_refine = 1
[]

[Variables]
  [./diffused]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
     type = SmoothCircleIC
     invalue = 1.0
     outvalue = 0.0
     radius = 10.0
     int_width = 2.0
     x1 = 25.0
     y1 = 25.0
    [../]
  [../]

  [./disp_x]
    order = FIRST
    family = LAGRANGE
  [../]

  [./disp_y]
    order = FIRST
    family = LAGRANGE
  [../]
[]

# Materials properties into AuxVariables - these are elemental variables, not nodal variables.
[AuxVariables]
  [./C11_aux]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./C12_aux]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./C13_aux]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./C14_aux]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./C15_aux]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./C16_aux]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./C22_aux]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./C23_aux]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./C24_aux]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./C25_aux]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./C26_aux]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./C33_aux]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./C34_aux]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./C35_aux]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./C36_aux]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./C44_aux]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./C45_aux]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./C46_aux]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./C55_aux]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./C56_aux]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./C66_aux]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Kernels]
  [./diff]
    type = Diffusion
    variable = diffused
  [../]
[]

[AuxKernels]
  [./matl_C11]
    type = MaterialSymmElasticityTensorAux
    tensor_matpro = elasticity_tensor
    index = 0
    variable = C11_aux
  [../]

  [./matl_C12]
    type = MaterialSymmElasticityTensorAux
    tensor_matpro = elasticity_tensor
    index = 1
    variable = C12_aux
  [../]

  [./matl_C13]
    type = MaterialSymmElasticityTensorAux
    tensor_matpro = elasticity_tensor
    index = 2
    variable = C13_aux
  [../]

  [./matl_C14]
    type = MaterialSymmElasticityTensorAux
    tensor_matpro = elasticity_tensor
    index = 3
    variable = C14_aux
  [../]

  [./matl_C15]
    type = MaterialSymmElasticityTensorAux
    tensor_matpro = elasticity_tensor
    index = 4
    variable = C15_aux
  [../]

  [./matl_C16]
    type = MaterialSymmElasticityTensorAux
    tensor_matpro = elasticity_tensor
    index = 5
    variable = C16_aux
  [../]

  [./matl_C22]
    type = MaterialSymmElasticityTensorAux
    tensor_matpro = elasticity_tensor
    index = 6
    variable = C22_aux
  [../]

  [./matl_C23]
    type = MaterialSymmElasticityTensorAux
    tensor_matpro = elasticity_tensor
    index = 7
    variable = C23_aux
  [../]

  [./matl_C24]
    type = MaterialSymmElasticityTensorAux
    tensor_matpro = elasticity_tensor
    index = 8
    variable = C24_aux
  [../]

  [./matl_C25]
    type = MaterialSymmElasticityTensorAux
    tensor_matpro = elasticity_tensor
    index = 9
    variable = C25_aux
  [../]

  [./matl_C26]
    type = MaterialSymmElasticityTensorAux
    tensor_matpro = elasticity_tensor
    index = 10
    variable = C26_aux
  [../]

  [./matl_C33]
    type = MaterialSymmElasticityTensorAux
    tensor_matpro = elasticity_tensor
    index = 11
    variable = C33_aux
  [../]

  [./matl_C34]
    type = MaterialSymmElasticityTensorAux
    tensor_matpro = elasticity_tensor
    index = 12
    variable = C34_aux
  [../]

  [./matl_C35]
    type = MaterialSymmElasticityTensorAux
    tensor_matpro = elasticity_tensor
    index = 13
    variable = C35_aux
  [../]

  [./matl_C36]
    type = MaterialSymmElasticityTensorAux
    tensor_matpro = elasticity_tensor
    index = 14
    variable = C36_aux
  [../]

  [./matl_C44]
    type = MaterialSymmElasticityTensorAux
    tensor_matpro = elasticity_tensor
    index = 15
    variable = C44_aux
  [../]

  [./matl_C45]
    type = MaterialSymmElasticityTensorAux
    tensor_matpro = elasticity_tensor
    index = 16
    variable = C45_aux
  [../]

  [./matl_C46]
    type = MaterialSymmElasticityTensorAux
    tensor_matpro = elasticity_tensor
    index = 17
    variable = C46_aux
  [../]

  [./matl_C55]
    type = MaterialSymmElasticityTensorAux
    tensor_matpro = elasticity_tensor
    index = 18
    variable = C55_aux
  [../]

  [./matl_C56]
    type = MaterialSymmElasticityTensorAux
    tensor_matpro = elasticity_tensor
    index = 19
    variable = C56_aux
  [../]

  [./matl_C66]
    type = MaterialSymmElasticityTensorAux
    tensor_matpro = elasticity_tensor
    index = 20
    variable = C66_aux
  [../]
[]


[Materials]
  [./test_material]
    type = LinearSingleCrystalPrecipitateMaterial
    block = 0
    disp_x = disp_x
    disp_y = disp_y
    C_matrix = '155.4 68.03 64.60 155.4 64.6 172.51 36.31 36.31 44.09'
    C_precipitate = '100 200 300 400 500 600  700 800 900'
    e_precipitate = '0.00551 0.0564 0.0570 0.0 0.0 0.0'
    n_variants = 1
    variable_names = 'diffused'
    all_21 = false
  [../]
[]

[BCs]
active = 'bottom top'
  [./bottom]
    type = DirichletBC
    variable = diffused
    boundary = '1'
    value = 1
  [../]

  [./top]
    type = DirichletBC
    variable = diffused
    boundary = '2'
    value = 0
  [../]
[]

[Executioner]
  type = Steady
  petsc_options = '-snes_mf_operator'
[]

[Output]
  file_base = Linear_Single_Xstal_Precip_Matl_out
  output_initial = true
  interval = 1
  exodus = true
  perf_log = true
[]
