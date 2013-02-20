[Mesh]
 #file = ScalingStudy_in.e
  type = GeneratedMesh
  dim = 3
  nx = 10
  ny = 10
  nz = 10
  xmin = 0
  xmax = 10
  ymin = 0
  ymax = 10
  zmin = 0
  zmax = 10
  elem_type = HEX8
[]

[Variables]
  [./concentration]
    order = THIRD
    family = HERMITE
      [./InitialCondition]
      type = ConstantIC
      value = 0.0562
    [../]
  [../]

  [./n1]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = SmoothCircleIC
      invalue = 1.6
      outvalue = 0
      radius = 1.8
      x1 = 2
      y1 = 2
      z1 = 2
    [../]
  [../]
  [./n2]
    order = FIRST
    family = LAGRANGE
  [../]
  [./n3]
    order = FIRST
    family = LAGRANGE
  [../]

  [./disp_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_z]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[AuxVariables]
  [./s11_aux]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./s12_aux]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./s13_aux]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./s22_aux]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./s23_aux]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./s33_aux]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./e11_aux]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./e12_aux]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./e13_aux]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./e22_aux]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./e23_aux]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./e33_aux]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./elem_ChemElastic_n1]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./elem_ChemElastic_n2]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./elem_ChemElastic_n3]
    order = CONSTANT
    family = MONOMIAL
  [../]

 [./elem_NucleationRate_n1]
    order = CONSTANT
    family = MONOMIAL
  [../]
 [./elem_NucleationRate_n2]
    order = CONSTANT
    family = MONOMIAL
  [../]
 [./elem_NucleationRate_n3]
    order = CONSTANT
    family = MONOMIAL
  [../]

 [./elem_NucleationProbability_n1]
    order = CONSTANT
    family = MONOMIAL
  [../]
 [./elem_NucleationProbability_n2]
    order = CONSTANT
    family = MONOMIAL
  [../]
 [./elem_NucleationProbability_n3]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[TensorMechanics]
  [./solid]
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
  [../]
[]

[Kernels]
  [./dcdt]
    type = TimeDerivative
    variable = concentration
  [../]

  [./dn1dt]
    type = TimeDerivative
    variable = n1
  [../]
  [./dn2dt]
    type = TimeDerivative
    variable = n2
  [../]
  [./dn3dt]
    type = TimeDerivative
    variable = n3
  [../]

  [./CHSolid]
    type = CHBulkPolyCoupled
    variable = concentration
    mob_name = M
    n_OP_variables = 3
    OP_variable_names = 'n1 n2 n3'
  [../]

  [./CHInterface]
    type = CHInterface
    variable = concentration
    kappa_name = kappa_c
    mob_name = M
    grad_mob_name = grad_M
  [../]

  [./ACSolidn1]
    type = ACBulkPolyCoupled
    variable = n1
    mob_name = L
    coupled_CH_var = concentration
    n_OP_vars = 3
    OP_var_names = 'n1 n2 n3'
    OP_number = 1
  [../]
  [./ACSolidn2]
    type = ACBulkPolyCoupled
    variable = n2
    mob_name = L
    coupled_CH_var = concentration
    n_OP_vars = 3
    OP_var_names = 'n1 n2 n3'
    OP_number = 2
  [../]
  [./ACSolidn3]
    type = ACBulkPolyCoupled
    variable = n3
    mob_name = L
    coupled_CH_var = concentration
    n_OP_vars = 3
    OP_var_names = 'n1 n2 n3'
    OP_number = 3
  [../]

  [./ACInterfacen1]
    type = ACInterface
    variable = n1
    mob_name = L
    kappa_name = kappa_n
  [../]
  [./ACInterfacen2]
    type = ACInterface
    variable = n2
    mob_name = L
    kappa_name = kappa_n
  [../]
  [./ACInterfacen3]
    type = ACInterface
    variable = n3
    mob_name = L
    kappa_name = kappa_n
  [../]

  [./ACTransformn1]
    type = ACTransformElasticDF
    variable = n1
    OP_number = 1
    OP_var_names = 'n1 n2 n3'
    n_OP_vars = 3
  [../]
  [./ACTransformn2]
    type = ACTransformElasticDF
    variable = n2
    OP_number = 2
    OP_var_names = 'n1 n2 n3'
    n_OP_vars = 3
  [../]
  [./ACTransformn3]
    type = ACTransformElasticDF
    variable = n3
    OP_number = 3
    OP_var_names = 'n1 n2 n3'
    n_OP_vars = 3
  [../]
[]

[AuxKernels]
  [./matl_s11]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 1
    index_j = 1
    variable = s11_aux
  [../]
 [./matl_s12]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 1
    index_j = 2
    variable = s12_aux
  [../]
  [./matl_s13]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 1
    index_j = 3
    variable = s13_aux
  [../]
  [./matl_s22]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 2
    index_j = 2
    variable = s22_aux
  [../]
  [./matl_s23]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 2
    index_j = 3
    variable = s23_aux
  [../]
  [./matl_s33]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 3
    index_j = 3
    variable = s33_aux
  [../]

  [./matl_e11]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 1
    index_j = 1
    variable = e11_aux
  [../]
 [./matl_e12]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 1
    index_j = 2
    variable = e12_aux
  [../]
  [./matl_e13]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 1
    index_j = 3
    variable = e13_aux
  [../]
  [./matl_e22]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 2
    index_j = 2
    variable = e22_aux
  [../]
  [./matl_e23]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 2
    index_j = 3
    variable = e23_aux
  [../]
  [./matl_e33]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 3
    index_j = 3
    variable = e33_aux
  [../]

  [./ChemElastic_n1]
    type = AuxChemElastic
    variable = elem_ChemElastic_n1
    coupled_conserved_var = concentration
    coupled_nonconserved_var = n1
    nonconserved_var_number = 1
    precip_conserved = 0.6
    precip_nonconserved = 1.6
  #  functional_c1 = 0.006
    execute_on = timestep
  [../]

  [./ChemElastic_n2]
    type = AuxChemElastic
    variable = elem_ChemElastic_n2
    coupled_conserved_var = concentration
    coupled_nonconserved_var = n2
    nonconserved_var_number = 2
    precip_conserved = 0.6
    precip_nonconserved = 1.6
  #  functional_c1 = 0.006
    execute_on = timestep
  [../]

 [./NucleationRate_n1]
    type = AuxNucleationRate
    variable = elem_NucleationRate_n1
    coupled_aux_var = elem_ChemElastic_n1
    #Kn1 = 0.008
    #Kn2 = 0.3

    gamma = 0.18
    scale_factor = 1

    Z = 0.1
    Beta_star = .0005
    linear_density = 5
    OP_var_names = 'n1 n2 n3'
    n_OP_vars = 3
    execute_on = timestep
  [../]
  [./NucleationRate_n2]
    type = AuxNucleationRate
    variable = elem_NucleationRate_n2
    coupled_aux_var = elem_ChemElastic_n2
    #Kn1 = 0.008
    #Kn2 = 0.3

    gamma = 0.18
    scale_factor = 1

    Z = 0.1
    Beta_star = .0005
    linear_density = 5
    OP_var_names = 'n1 n2 n3'
    n_OP_vars = 3
    execute_on = timestep
  [../]
  [./NucleationRate_n3]
    type = AuxNucleationRate
    variable = elem_NucleationRate_n3
    coupled_aux_var = elem_ChemElastic_n3
    #Kn1 = 0.008
    #Kn2 = 0.3

    gamma = 0.18
    scale_factor = 1

    Z = 0.1
    Beta_star = .0005
    linear_density = 5
    OP_var_names = 'n1 n2 n3'
    n_OP_vars = 3
    execute_on = timestep
  [../]

 [./AuxNucleationProbability_n1]
    type = AuxNucleationProbability
    variable = elem_NucleationProbability_n1
    coupled_aux_var = elem_NucleationRate_n1
    coupled_variables = 'n1 n2 n3'
    n_OP_vars = 3
    execute_on = timestep
 [../]
 [./AuxNucleationProbability_n2]
    type = AuxNucleationProbability
    variable = elem_NucleationProbability_n2
    coupled_aux_var = elem_NucleationRate_n2
    coupled_variables = 'n1 n2 n3'
    n_OP_vars = 3
    execute_on = timestep
 [../]
 [./NucleationProbability_n3]
    type = AuxNucleationProbability
    variable = elem_NucleationProbability_n3
    coupled_aux_var = elem_NucleationRate_n3
    coupled_variables = 'n1 n2 n3'
    n_OP_vars = 3
    execute_on = timestep
  [../]
[]

[BCs]
  [./disp_x_BC]
    type = DirichletBC
    variable = disp_x
    boundary = '0 1 2 3 4 5'
    value = 0.0
  [../]

  [./disp_y_BC]
    type = DirichletBC
    variable = disp_y
    boundary = '0 1 2 3 4 5'
    value = 0.0
  [../]

  [./disp_z_BC]
    type = DirichletBC
    variable = disp_z
    boundary = '0 1 2 3 4 5'
    value = 0.0
  [../]

#  [./c_BC]
#    type = NeumannBC
#    variable = concentration
#    boundary = '0 1 2 3'
#    value = 0.0
#  [../]

 # [./n1_BC]
 #   type = NeumannBC
 #   variable = n1
 #   boundary = '0 1 2 3'
 #   value = 0.0
 # [../]

 # [./n2_BC]
 #   type = NeumannBC
 #   variable = n2
 #   boundary = '0 1 2 3'
 #   value = 0.0
 # [../]

 # [./n3_BC]
 #   type = NeumannBC
 #   variable = n3
 #   boundary = '0 1 2 3'
 #   value = 0.0
 # [../]
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
    disp_z = disp_z
    #reading C_11  C_12  C_13  C_22  C_23  C_33  C_44  C_55  C_66
    C_ijkl ='155.4 68.03 64.6 155.4  64.6 172.51 36.31 36.31 44.09'
    C_precipitate ='155.4 68.03 64.6 155.4  64.6 172.51 36.31 36.31 44.09'
    #reading        S_11   S_22  S_33 S_23 S_13 S_12
    e_precipitate = '0.00551  0.0564  0.0570  0.0  0.0  0.0'
    n_variants = 3
    variable_names = 'n1 n2 n3'
    all_21 = false
  [../]
[]

[Postprocessors]
  [./ElementIntegral_n1]
    output = file
    type = ElementIntegralVariablePostprocessor
    variable = n1
  [../]

  [./ElementIntegral_n2]
    output = file
    type = ElementIntegralVariablePostprocessor
    variable = n2
  [../]

  [./ElementIntegral_n3]
    output = file
    type = ElementIntegralVariablePostprocessor
    variable = n3
  [../]

  [./ElementIntegral_c]
    output = file
    type = ElementIntegralVariablePostprocessor
    variable = concentration
  [../]

  [./NodalMaxValue_n1]
    output = file
    type = NodalMaxValue
    variable = n1
  [../]

  [./NodalMaxValue_n2]
    output = file
    type = NodalMaxValue
    variable = n2
  [../]

  [./NodalMaxValue_n3]
    output = file
    type = NodalMaxValue
    variable = n3
  [../]

  [./NodalMaxValue_c]
    output = file
    type = NodalMaxValue
    variable = concentration
  [../]

  [./ndofs]
    type = PrintDOFs
    output = file
  [../]
[]

[UserObjects]
  [./NLUO]
    type = NucleationLocationUserObject
    coupled_aux_vars = 'elem_NucleationProbability_n1 elem_NucleationProbability_n2 elem_NucleationProbability_n3'
    n_coupled_aux = 3
    dwell_time = 0.1
    num_orientations = 3
    execute_on = timestep
  [../]

  [./NISM]
    type = NucleusIntroductionSolutionModifier
    nucleation_userobject = NLUO
    variables = 'n1 n2 n3'
    num_vars = 3
    seed_value = 1.6
    radius = 1.8
    int_width = 0.9
    execute_on = custom
  [../]
[]

[Executioner]
  type = MeshSolutionModify
  scheme = 'crank-nicolson'

  num_steps = 10
  dt = 0.01
  dtmin = 0.0001
  dtmax = 0.1
  percent_change = 0.1

  start_time = 0.0
 # end_time = 100

  abort_on_solve_fail = true
  adapt_nucleus = 5
  adapt_cycles = 1

  use_nucleation_userobject = true
  nucleation_userobject = NLUO

  petsc_options = -snes_mf_operator
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]

#[Adaptivity]
#   marker = combo
# [./Markers]
#    [./NM]
#      type = NucleationMarker
#      nucleation_userobject = NLUO
#      max_h_level = 5
#    [../]
#    [./EFMHM]
#      type = ErrorFractionMaxHMarker
#      coarsen = 0.05
#      refine = 0.75
#      indicator = GJI
#      max_h_level = 5
#    [../]
#    [./combo]
#      type = ComboMarker
#      markers = 'NM EFMHM'
#    [../]
# [../]

# [./Indicators]
#   [./GJI]
#    type = GradientJumpIndicator
#    variable = concentration
#   [../]
# [../]
#[]

[Output]
  file_base = ScalingStudy
  output_initial = true
  interval = 1
  exodus = true
  perf_log = true
  postprocessor_csv = true
[]