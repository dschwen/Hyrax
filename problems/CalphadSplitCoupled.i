[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 50 #100
  ny = 50 #100
  nz = 50 #100
  xmin = 0
  xmax = 50E-9
  ymin = 0
  ymax = 50E-9
  zmin = 0
  zmax = 50E-9
  elem_type = QUAD4
  #elem_type = QUAD9
  #elem_type = HEX8
[]

[Variables]
  [./concentration]
   # order = SECOND
    order = FIRST
    family = LAGRANGE
  [../]

  [./mu]
    order = FIRST
    #order = SECOND
    family = LAGRANGE
  [../]

  [./n]
    order = FIRST
    family = LAGRANGE
  [../]

  [./temperature]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[AuxVariables]
  [./Galpha]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./Gdelta]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./dGalphadc]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./dGdeltadc]
    order = CONSTANT
    family = MONOMIAL
  [../]


  [./d2Galphadc2]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./d2Gdeltadc2]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./auxGalpha]
    type = MaterialRealAux
    variable = Galpha
    property = G_AB1CD1
  [../]

  [./auxGdelta]
    type = MaterialRealAux
    variable = Gdelta
    property = G_AB1CD2
  [../]

  [./auxdGalphaDc]
    type = MaterialRealAux
    variable = dGalphadc
    property = dGAB1CD1_dc
  [../]

  [./auxdGdeltaDc]
    type = MaterialRealAux
    variable = dGdeltadc
    property = dGAB1CD2_dc
  [../]

  [./auxd2Galphadc2]
    type = MaterialRealAux
    variable = d2Galphadc2
    property = d2GAB1CD1_dc2
  [../]

   [./auxd2Gdeltadc2]
    type = MaterialRealAux
    variable = d2Gdeltadc2
    property = d2GAB1CD2_dc2
   [../]
[]

[ICs]
  [./PSSCIC]
      #type = SmoothCircleIC
      #type = ConstantIC
      type = RandomIC
      #value = 0.66
      variable = concentration
      min = 0.6549
      max = 0.6551
      #int_width = 3E-9
      #invalue = 0.6
      #outvalue = 0.1
      #radius = 10E-9
      #x1 = 25E-9
      #y1 = 25E-9
  [../]

  [./PSSCIC_n]
      #type = SmoothCircleIC
      variable = n
      #int_width = 3E-9
      #invalue = 1
      #outvalue = 0
      #radius = 10E-9
      #x1 = 25E-9
      #y1 = 25E-9
      type = RandomIC
      min = 0.0
      max = 0.001
  [../]

  [./t_init]
    type = ConstantIC
    variable = temperature
    value = 600
  [../]
[]

[Preconditioning]
#active = 'SMP'
#  [./PBP]
#   type = PBP
#   solve_order = 'w c'
#   preconditioner = 'AMG ASM'
#   off_diag_row = 'c '
#   off_diag_column = 'w '
#  [../]

  [./SMP]
   type = SMP
   off_diag_row = 'mu concentration'
   off_diag_column = 'concentration mu'
  [../]
[]

[Kernels]
  [./dcdt]
    type = CoupledImplicitEuler
    variable = mu
    v = concentration
  [../]

  [./mu_residual]
    type = SplitCHWRes
    variable = mu
    mob_name = M
    c = concentration
  [../]

  [./conc_residual]
    type = CHCoupledCalphadSplit
    variable = concentration
    kappa_name = kappa_c
    coupled_OP_var = n
    w = mu
    #Well_height = 1
  [../]

  [./dndt]
    type = TimeDerivative
    variable = n
  [../]

  [./ACSolidn]
    type = ACCoupledCalphad
    variable = n
    mob_name = L
    #coupled_CH_var = concentration
    n_OP_vars = 1
    OP_var_names = 'n'
    OP_number = 1
  [../]

  [./ACInterfacen1]
    type = ACInterface
    variable = n
    mob_name = L
    kappa_name = kappa_n
  [../]

  [./dTdt]
    type = TimeDerivative
    variable = temperature
  [../]


  [./THeat]
    type = Heat
    variable = temperature
  [../]
[]


[Materials]
  [./calphad]
    type = ZrHCalphad
    block = 0

    #H_Zr_D0 = 7.00e5 #um^2/s
    #H_ZrH2_D0 = 1.53e5	 # um^2/s
    #H_Zr_Q0 =  4.456e4 #J/mol
    #H_ZrH2_Q0 = 5.885E4 #J/mol

    mobility_AC = 2E7
    mobility_CH = 2E-10

    kappa_CH = 7.812E-10
    kappa_AC = 7.812E-9

    well_height = 60000

    molar_volume = 1.4E-5

    thermal_diffusivity = 1

    coupled_temperature = temperature
  [../]

  [./alphaZr]
   type = CalphadAB1CD1
   block = 0

   pure_endpoint_low_coeffs = '-7827.595
                                 125.64905
                                 -24.1618
                                  -0.00437791
                               34971.0' #HCP_Zr
  #this is fucked...

 pure_endpoint_high_coeffs = '2589675
                                 -34719.21
                                   5126.9713
                                     -3.250187
                             -208070050'  #H2_gas
   mixture_coeffs = '-45965
                         41.6
                          0'  #FCC_ZrH
   L0_coeffs = '0 0'
   L1_coeffs = '0 0'
   low_coeffs = '   -26549.141703672653
                 -31205591.610163368
                       295.24365293254311
                7516799743.7832689
                        -0.42527422354460137
                     -4439.0959906994285
                         0.00022153720566517584
                    -75983.150207815517
                  27656668.097620554
                       -18.137337020144734'

   high_coeffs = '-1.1748697700305443E7
                   7.1428722799052447E7
                   2.0711053625846262E2
                  -1.4437278547848886E8
                   3.5284806111023137E-1
                   9.7186016495869264E7
                  -1.9493853217896482E-4
                  -1.9533177347714059E3
                   2.0853116312248094E3
                  -3.9319750123563216E-3'

    coupled_temperature = temperature
    coupled_concentration = concentration
  [../]

  [./deltaZrH2]
   type = CalphadAB1CD2
   block = 0
   pure_endpoint_low_coeffs = '-227.595
                               124.74905
                               -24.1618
                                -0.00437791
                             34971' #FCC_Zr
    pure_endpoint_high_coeffs = '2589675
                                 -34719.21
                                   5126.9713
                                     -3.250187
                             -208070050'  #H2_gas
   mixture_coeffs =  '-170490
                          208.2
                           -9.47' #FCC_ZrH2'
   L0_coeffs = '14385 -6.0'
   L1_coeffs = '-106445 87.3'

   low_coeffs =  '2.6978175664339593E4
                 -5.2267355124714538E7
                  7.7077694451557377E1
                  2.5671915976342548E10
                 -3.0872406268790027E-2
                 -4.5394250810591288E1
                 -1.5479208178631023E-5
                 -8.3176435859210935E4
                  4.4517036341637028E6
                  2.1739434332085402E1'

    high_coeffs = ' -4.8353171133552969E7
                     2.2151622968687481E8
                    -6.0924808041164090E2
                    -3.3818504692030811E8
                     5.3781670052188990E-1
                     1.7206754571759129E8
                    -2.6254489645679825E-4
                     7.6102442020895887E2
                    -4.4943040350977139E2
                    -9.0858040404436430E-2'

    coupled_temperature = temperature
    coupled_concentration = concentration

    pure_EP1_phase1_coeffs = '-7827.595
                                 125.64905
                                 -24.1618
                                  -0.00437791
                               34971.0' #HCP_Zr
  [../]
[]

[Postprocessors]
  [./ElementIntegral_c]
    type = ElementIntegralVariablePostprocessor
    variable = concentration
    output = file
  [../]

  [./ElementIntegral_n]
    type = ElementIntegralVariablePostprocessor
    variable = n
    output = file
  [../]

  [./max_c]
    type = NodalMaxValue
    variable = concentration
    output = file
  [../]

  [./max_n]
    type = NodalMaxValue
    variable = n
    output = file
  [../]

  [./volume_fraction]
    type = NodalVolumeFraction
    variable = 'n'
    use_single_map = true
    threshold = 0.75
    execute_on = timestep
    mesh_volume = Volume
  [../]

  [./Volume]
    type = VolumePostprocessor
    execute_on = initial
  [../]
[]

[Executioner]
  type = Transient
  #scheme = 'crank-nicolson'
  scheme = 'BDF2'

  #Preconditioned JFNK (default)
  #solve_type = 'PJFNK'
  solve_type = 'JFNK'
  #solve_type = 'NEWTON'

  #petsc_options_iname = '-pc_type'
  #petsc_options_value = 'lu'

  l_max_its = 20
  #l_tol = 1.0e-5

  #nl_max_its = 40
  #nl_rel_tol = 1.0e-8
  nl_max_its = 20
  nl_abs_tol = 1.25E-4

  start_time = 0.0
  num_steps = 2
  dt = 1E-20
  dtmin = 1E-22
[]

[Output]
  file_base = CHCoupledCalphadSplit
  output_initial = true
  interval = 1
  linear_residuals = true
  exodus = true
  perf_log = true
  all_var_norms = true

  postprocessor_csv = true
[]
