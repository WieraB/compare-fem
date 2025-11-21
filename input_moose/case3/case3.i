#-------------------------------------------------------------------------
# pyvale: gmsh,3Dstc,1mat,thermal,steady,
#-------------------------------------------------------------------------

#-------------------------------------------------------------------------
#_* MOOSEHERDER VARIABLES - START

# Thermal Loads/BCs
coolantTemp = 100.0      # degC
heatTransCoeff = 125.0e3 # W.m^-2.K^-1
surfHeatFlux = 5.0e6    # W.m^-2

# Material Properties: Pure (OFHC) Copper at 250degC
cuDensity = 8829.0  # kg.m^-3
cuThermCond = 384.0 # W.m^-1.K^-1
cuSpecHeat = 406.0  # J.kg^-1.K^-1

#** MOOSEHERDER VARIABLES - END
#-------------------------------------------------------------------------

[Mesh]
    type = FileMesh
    file = 'mesh_block_pipe_refined.msh'
[]

[Variables]
    [temperature]
        initial_condition = ${coolantTemp}
    []
[]

[Kernels]
    [heat_conduction]
        type = HeatConduction
        variable = temperature
    []
[]

[Materials]
    [copper_thermal]
        type = HeatConductionMaterial
        thermal_conductivity = ${cuThermCond}
        specific_heat = ${cuSpecHeat}
    []
    [copper_density]
        type = GenericConstantMaterial
        prop_names = 'density'
        prop_values = ${cuDensity}
    []
[]

[BCs]
    [heat_flux_out]
        type = ConvectiveHeatFluxBC
        variable = temperature
        boundary = 'bc-pipe-htc'
        T_infinity = ${coolantTemp}
        heat_transfer_coefficient = ${heatTransCoeff}
    []
    [heat_flux_in]
        type = NeumannBC
        variable = temperature
        boundary = 'bc-top-heatflux'
        value = ${surfHeatFlux}
    []
[]

[Preconditioning]
    [SMP]
        type = SMP
        full = true
    []
[]

[Executioner]
    type = Steady
    solve_type = 'NEWTON'
    l_max_its = 100       # default = 1000
    nl_max_its = 100
    l_tol = 1e-5          # default = 1e-5, worse match with ngsove if it's 1e-6
    nl_abs_tol = 1e-6     # default = 1e-50, set 1e-6
    nl_rel_tol = 1e-6     # default = 1e-8, set 1e-6
    petsc_options_iname = '-pc_type -pc_hypre_type'
    petsc_options_value = 'hypre boomeramg'
[]

[Outputs]
    exodus = true
[]