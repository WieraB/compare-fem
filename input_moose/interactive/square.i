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

[Executioner]
    type = Steady
[]

[Outputs]
    exodus = true
[]