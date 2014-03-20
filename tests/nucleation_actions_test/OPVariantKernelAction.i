# this input file is to test the Actions created for concurrent nucleation and growth with 3 orientation variants

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 40
  ny = 40
  nz = 0
  xmin = 0
  xmax = 6
  ymin = 0
  ymax = 6
  zmin = 0
  zmax = 0
  elem_type = QUAD4
[]

[Variables]
  [./concentration]
    order = THIRD
    family = HERMITE
    [./InitialCondition]
      type = SmoothCircleIC
      radius = 1
      int_width = 0.9
      x1 = 3
      y1 = 3
      invalue = 0.6
      outvalue = 0.05
    [../]
  [../]

  [./n1]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = SmoothCircleIC
      radius = 1
      int_width = 0.9
      x1 = 3
      y1 = 3
      invalue = 1.6
      outvalue = 0.0
    [../]
  [../]

  [./n2]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 0.001
      max = 0.01
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

[Kernels]
  [./dcdt]
    type = TimeDerivative
    variable = concentration
  [../]

  [./CHSolid]
   type = CHBulkPolyCoupled
    variable = concentration
    mob_name = M
    n_OP_variables = 2
    OP_variable_names = 'n1 n2'
  [../]

  [./CHInterface]
    type = CHInterface
    variable = concentration
    kappa_name = kappa_c
    mob_name = M
    grad_mob_name = grad_M
  [../]
[]

[TensorMechanics]
  [./solid]
    disp_x = disp_x
    disp_y = disp_y
  [../]
[]

[OPVariantKernel]
    number_OPs = 2
    OP_name_base = n

    #for ACBulkPolyCoupled
    coupled_CH_var = concentration

    #for ACInterface
    kappa_name_OP = kappa_n

    #for TimeDerivative
[]

[BCs]
  [./c_BC]
    type = NeumannBC
    variable = concentration
    boundary = '1'
    value = 1
  [../]

  [./fixed_x_BC]
    type = DirichletBC
    variable = disp_x
    boundary = '0 1 2 3'
    value = 0.0
  [../]
  [./fixed_y_BC]
    type = DirichletBC
    variable = disp_y
    boundary = '0 1 2 3'
    value = 0.0
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

  [./test_material]
    type = LinearSingleCrystalPrecipitateMaterial
    block = 0
    disp_x = disp_x
    disp_y = disp_y
    #reading C_11  C_12  C_13  C_22  C_23  C_33  C_44  C_55  C_66
    C_ijkl ='155.4 68.03 64.6 155.4  64.6 172.51 36.31 36.31 44.09'
    C_precipitate ='155.4 68.03 64.6 155.4  64.6 172.51 36.31 36.31 44.09'
    #reading        S_11   S_22  S_33 S_23 S_13 S_12
    e_precipitate = '0.00551  0.0564  0.0570  0.0  0.0  0.0'
    n_variants = 2
    variable_names = 'n1 n2'
    all_21 = false
  [../]
[]

[Executioner]
  type = Transient
  scheme = 'bdf2'

  [./TimeStepper]
    type = ConstantDT
    dt = 0.01
  [../]

  num_steps = 1

  l_max_its = 50

  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'

  petsc_options = '' #-ksp_monitor'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]

[Outputs]
  file_base = OPVariantKernelAction
  output_initial = true
  exodus = true
  [./console]
    type = Console
    perf_log = true
    linear_residuals = true
  [../]
[]
