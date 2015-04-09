# this input file is to test the concurrent nucleation and growth with mesh adaptivity

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 10
  ny = 10
  nz = 0
  xmin = 0
  xmax = 76.8 #0.3*256
  ymin = 0
  ymax = 76.8
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
      min = 0.0
      max = 0.1
    [../]
  [../]

  [./n]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 0.0
      max = 1E-5
    [../]
  [../]
[]

[AuxVariables]
  [./elemental_Supersaturation]
    order = CONSTANT
    family = MONOMIAL
   [../]

 [./elemental_NucleationProbability]
    order = CONSTANT
    family = MONOMIAL
  [../]

 [./elemental_NucleationRate]
    order = CONSTANT
    family = MONOMIAL
  [../]

 [./temperature]
   order = CONSTANT
   family = MONOMIAL
 [../]
[]

[Kernels]
  [./dcdt]
    type = TimeDerivative
    variable = concentration
  [../]

  [./diff_c]
    type = Diffusion
    variable = concentration
  [../]

 [./dndt]
    type = TimeDerivative
    variable = n
  [../]

  [./diff_n]
    type = Diffusion
    variable = n
 [../]

#  [./dn1dt]
#    type = TimeDerivative
#    variable = n1
#  [../]

#  [./CHSolid]
#    type = CHBulkCoupled
#    variable = concentration
#    mob_name = M
#    coupled_OP_var = n1
#  [../]

#  [./CHInterface]
#    type = CHInterface
#    variable = concentration
#    kappa_name = kappa_c
#    mob_name = M
#    grad_mob_name = grad_M
#  [../]

#  [./ACSolidn1]
#    type = ACBulkCoupled
#    variable = n1
#    mob_name = L
#    coupled_CH_var = concentration
#  [../]

#  [./ACInterfacen1]
#    type = ACInterface
#    variable = n1
#    mob_name = L
#    kappa_name = kappa_n
#  [../]
[]

[AuxKernels]
  [./Supersaturation]
    type = AuxSupersaturation
    variable = elemental_Supersaturation
    coupled_var = concentration
    functional_c1 = 0.006
    execute_on = timestep_end
  [../]

  [./NucleationRate]
    type = AuxVolumetricFullNucleationRate
    variable = elemental_NucleationRate
    coupled_bulk_energy_change = elemental_Supersaturation
    OP_variable_names = 'n'
    n_OP_variables = 1
    T = temperature
    X = concentration

    gamma = 1
    Kb = 1
    linear_density = 1
    jump_distance = 1
    rate_volume = 1

    execute_on = timestep_end
    length_scale_factor = 1
    scale_factor = 5E14
  [../]

  [./NucleationProbability]
    type = AuxNucleationProbability
    variable = elemental_NucleationProbability
    coupled_aux_var = elemental_NucleationRate
    coupled_variables = 'n'
    n_OP_vars = 1
    OP_threshold = 0.0001
    execute_on = timestep_end
  [../]

  [./auxtemp]
    type = AuxTemperature
    variable = temperature
    temp_in_K = 600
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

    #interface energies are scaled
    kappa_CH = 1.7708
    kappa_AC = 17.708

    #well height and molar volume remain unscaled.
    well_height = 40500
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
                          0'  #FCC_ZrH
   L0_coeffs = '0 0'
   L1_coeffs = '0 0'


    coupled_temperature = temperature
    coupled_concentration = concentration
  [../]

  [./deltaZrH2]
   type = CalphadAB1CD2Material
   block = 0
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

[BCs]
#  [./c_BC]
#    type = NeumannBC
#    variable = concentration
#    boundary = '0 1 2 3'
#    value = 0.0
#  [../]

#  [./n1_BC]
#    type = NeumannBC
#    variable = n1
#    boundary = '0 1 2 3'
#    value = 0.0
#  [../]
[]

#[Materials]
#  [./constant]
#    type = PFMobilityLandau
#    block = 0
#    mob_CH = 0.4
#    mob_AC = 0.4
#    kappa_CH = 1.5
#    kappa_AC = 1.5
#    A1 = 18.5
#    A2 = -8.5
#    A3 = 11.5
#    A4 = 4.5
#    C1 = 0.006
#    C2 = 0.59
#  [../]
#[]

[UserObjects]
  [./NLUO]
    type = NucleationLocationUserObject
    #variable = concentration
    coupled_aux_vars = 'elemental_NucleationProbability'
    n_coupled_aux = 1
    dwell_time = 0.1
    num_orientations = 1
    boundary_width = 5
    execute_on = timestep_end
  [../]

  [./NISM]
    type = NucleusIntroductionSolutionModifier
    nucleation_userobject = NLUO
    variables = 'n'
    num_vars = 1
    seed_value = 1.6
    radius = 1.8
    int_width = 0.9
    execute_on = custom
  [../]
[]

[Executioner]
  type = MeshSolutionModify
  scheme = 'crank-nicolson'
  num_steps = 3
  [./TimeStepper]
    type = ConstantDT #SolutionTimeAdaptiveDT
    dt = 0.1
  [../]

  adapt_nucleus = 4


  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'

  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
  use_nucleation_userobject = true
  nucleation_userobject = NLUO
[]

[Adaptivity]
   marker = NM
 [./Markers]
    [./NM]
      type = NucleationMarker
      nucleation_userobject = NLUO
      max_h_level = 4
    [../]
 [../]
[]

[Outputs]
  file_base = AMR_volumetric_rate_test
  exodus = true
  output_initial = true
  [./console]
    type = Console
    perf_log = true
    linear_residuals = true
  [../]
  checkpoint = false
[../]
