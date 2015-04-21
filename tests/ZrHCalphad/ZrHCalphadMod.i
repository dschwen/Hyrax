#This simulation is non-dimensionalized

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 2 
  ny = 2 
  nz = 0 
  xmin = 0 
  xmax = 2 
  ymin = 0 
  ymax = 2 
  zmin = 0
  zmax = 0
  elem_type = QUAD4
[]

[Variables]
  [./concentration]
    order = FIRST
    family = LAGRANGE
  [../]

  [./mu]
    order = FIRST
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
  [./CIC_c]
      type = ConstantIC
      variable = concentration
      value = 0.1
  [../]

  [./RIC_n]
      type = RandomIC
      variable = n
      min = 0
      max = 0.001
  [../]

  [./t_init]
    type = ConstantIC
    variable = temperature
    value = 800
  [../]
[]

[Preconditioning]
 active = 'SMP'
  [./SMP]
   type = SMP
   off_diag_row = 'concentration concentration concentration mu mu mu n n n temperature temperature temperature'
   off_diag_column = 'mu n temperature concentration n temperature concentration mu temperature concentration mu n'
  [../]
[]

[Kernels]
  [./dcdt]
    type = CoupledImplicitEuler
    variable = mu
    v = concentration
  [../]

  [./mu_residual]
    type = SplitCoupledCHWRes
    variable = mu
    mob_name = M
    c = concentration
    T = temperature
    n_OP_vars = 1
    OP_var_names = 'n'
  [../]

  [./conc_residual]
    type = CHCoupledCalphadSplit
    variable = concentration
    kappa_name = kappa_c
    w = mu
    n_OP_vars = 1
    OP_var_names = 'n'
    T = temperature
    #Well_height = 1
    scaling_factor = 2.49410145E-9
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
    w = mu
    T = temperature
    c = concentration
    scaling_factor = 2.49410145E-9
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
    w = mu
    c = concentration
    OP_var_names = 'n'
    n_OP_vars = 1
  [../]
[]


[Materials]
  [./calphad]
    type = ZrHCalphadDiffusivity
    block = 0

    OP_variable = n
    concentration = concentration

    H_Zr_D0 = 7.00e-7    #m^2/s
    H_ZrH2_D0 = 1.53e-7  # m^2/s
    H_Zr_Q0 =  4.456e4   #J/mol
    H_ZrH2_Q0 = 5.885E4  #J/mol

   #ALWAYS have mobility_AC = 1 for nondimensionalization (L*)!!
    #we pick L = 1E0, and so L* = 1
    mobility_AC = 1E0
    #to scale mobility, M* = M/L*l^2
    CH_mobility_scaling = 1E-18

    kappa_CH = 1.9484
    kappa_AC = 19.484

    #well height and molar volume remain unscaled.
    well_height = 0

    molar_volume = 1.4E-5

    thermal_diffusivity = 1

    coupled_temperature = temperature
  [../]

  [./alphaZr]
   type = CalphadAB1CD1Material
   block = 0

   pure_endpoint_low_coeffs = '-7827.595
                                 125.64905
                                 -24.1618
                                  -0.00437791
                               34971.0' #HCP_Zr

  pure_endpoint_high_coeffs = '8055.336
                                 -243.791
                                   18.3135
                                   -0.034513
                              -734182.8'  #H2_gas

   mixture_coeffs = '-45965
                         41.6
                          0'  #HCP_ZrH
   L0_coeffs = '0 0'
   L1_coeffs = '0 0'

   coupled_temperature = temperature
   coupled_concentration = concentration
  [../]

  [./deltaZrH2]
   type = CalphadAB1CD2ModMaterial
   block = 0

   mod_cutoff = 0.57
   A1 = '18.27 -12934'
   A2 = '31.73 -26537'

   A1_to_0 = 707.936
   A2_to_0 = 836

   pure_endpoint_low_coeffs = '-227.595
                               124.74905
                               -24.1618
                                -0.00437791
                             34971' #FCC_Zr

  pure_endpoint_high_coeffs = '8055.336
                                 -243.791
                                   18.3135
                                   -0.034513
                              -734182.8'  #H2_gas
   mixture_coeffs =  '-170490
                          208.2
                           -9.47' #FCC_ZrH2'

   L0_coeffs = '14385 -6.0'
   L1_coeffs = '-106445 87.3'


    coupled_temperature = temperature
    coupled_concentration = concentration

    pure_EP1_phase1_coeffs = '-7827.595
                                 125.64905
                                 -24.1618
                                  -0.00437791
                               34971.0' #HCP_Zr
  [../]
[]

[Executioner]
  type = Transient
  scheme = 'BDF2'

  [./TimeStepper]
    type = SolutionTimeAdaptiveDT
    dt = 5e-5
    percent_change = 0.1
  [../]

  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'

  petsc_options_iname = '-pc_type -ksp_gmres_restart'
  petsc_options_value = 'ilu 50'

   l_max_its = 100
   #l_tol = 1.0e-5

  #nl_max_its = 40
  nl_rel_tol = 1e-7
  #nl_max_its = 20
  #nl_abs_tol = 2.15E-6

  start_time = 0.0
  num_steps = 2
  dtmax = 1E-1
  dtmin = 1E-7
[]

[Outputs]
  file_base = ZrHCalphadMod
  output_initial = true
  exodus = true
  [./console]
    type = Console
    perf_log = true
    linear_residuals = true
  [../]
  csv = true
[]
