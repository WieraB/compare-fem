import time
from pathlib import Path
from typing import Any
import dataclasses
import shutil
import numpy as np
from scipy.spatial import KDTree
import pyvista as pv

#pyvale imports
import pyvale.dataset as dataset
import pyvale.sensorsim as sens
from pyvale.mooseherder import (MooseConfig,
                                MooseRunner,
                                ExodusReader)
#%%
# Set file locations
case_name = "case3"
mesh_name = "mesh_block_pipe_refined.msh"

mesh_path = "./meshes/"
output_path = f"./output/{case_name}/{case_name}_moose.vtu"
input_file_path = f"./input_moose/{case_name}/{case_name}.i"
input_path = f"./input_moose/{case_name}/"

mesh_file_path = mesh_path + mesh_name
mesh_file_path_moose = input_path + mesh_name
shutil.copyfile(mesh_file_path, mesh_file_path_moose)

config = {'main_path': Path.home()/ 'moose',
          'app_path': Path.home() / 'proteus',
          'app_name': 'proteus-opt'}
moose_config = MooseConfig(config)

moose_runner = MooseRunner(moose_config)

moose_runner.set_run_opts(n_tasks = 1,
                          n_threads = 1,
                          redirect_out = False)

moose_input = Path(input_file_path) # It actually works now :)

moose_runner.set_input_file(moose_input)

print(moose_runner.get_arg_list())
print()

#%%
# Run the simulation

start_time = time.perf_counter()
moose_runner.run()
run_time = time.perf_counter() - start_time

print()
print("-"*80)
print(f'MOOSE run time = {run_time:.3f} seconds')
print("-"*80)
print()

#%%
# Save the results
# WARNING!!! If there is an error, it is likely because the version of .msh file is above 2 for pyvista.

input_path = Path(input_path)

output_exodus = input_path / (moose_input.stem + "_out.e")
exodus_reader = ExodusReader(output_exodus)

print("\nReading exodus file with ExodusReader:")
print(output_exodus.resolve())
print()

read_config = exodus_reader.get_read_config()
sens.SimTools.print_dataclass_fields(read_config)

read_config.time = False
read_config.coords = True
read_config.connect = False
sim_data = exodus_reader.read_sim_data(read_config)

sens.SimTools.print_sim_data(sim_data)

sol_unordered = sim_data.node_vars['temperature'][:, 1]
points_unordered = sim_data.coords

res = pv.read(mesh_file_path)
points = res.points

tree2 = KDTree(points_unordered)

_, indices = tree2.query(points)

sol = sol_unordered[indices]

res["sol"] = sol

res.save(output_path)