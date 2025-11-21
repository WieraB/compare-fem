#-------------------------------------------------------------------------
# pyvale: simple,2Dplate,1mat,thermal,steady,
#-------------------------------------------------------------------------

#-------------------------------------------------------------------------
#_* MOOSEHERDER VARIABLES - START

# Thermal Loads/BCs
coolantTemp = 0.0      # degC

# Material Properties
cuDensity = 1.0  # kg.m^-3
cuThermCond = 1.0 # W.m^-1.K^-1
cuSpecHeat = 1.0  # J.kg^-1.K^-1

#** MOOSEHERDER VARIABLES - END
#-------------------------------------------------------------------------

[Mesh]
    type = FileMesh
    file = 'mesh_square.msh'
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
    [bc_left]
        type = DirichletBC
        variable = temperature
        boundary = 'Left'
        value = '0'
    []
    [heat_right]
        type = DirichletBC
        variable = temperature
        boundary = 'Right'
        value = '0'
    []
    [bc_top]
        type = FunctionNeumannBC
        variable = temperature
        boundary = 'Top'
        function = '10*sin(6*x)'
    []
    [bc_bottom]
        type = FunctionNeumannBC
        variable = temperature
        boundary = 'Bottom'
        function = '10*sin(6*x)'
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