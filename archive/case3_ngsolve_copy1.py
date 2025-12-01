from ngsolve import *
from mpi4py import MPI
from netgen.meshing import Mesh
from netgen.occ import unit_square
import meshio
from ngsolve import *
from ngsolve.webgui import Draw
import netgen
from netgen.geom2d import SplineGeometry
from netgen.read_gmsh import ReadGmsh
import time
import numpy as np
import pyvista as pv

#%%
# Set file locations
case_name = "case3"
mesh_path = "./meshes/mesh_block_pipe_refined.msh"

output_path = f"./output/{case_name}/{case_name}_ngsolve.vtu"

#%%
# Define material properties and loads

# Thermal Loads/BCs
coolantTemp = 100.0      # degC
heatTransCoeff = 125.0e3 # W.m^-2.K^-1
surfHeatFlux = 5.0e6    # W.m^-2
T_inf = coolantTemp

# Material Properties: Pure (OFHC) Copper at 250degC
cuDensity = 8829.0  # kg.m^-3
cuThermCond = 384.0 # W.m^-1.K^-1
cuSpecHeat = 406.0  # J.kg^-1.K^-1

#%%
# Load the mesh
# WARNING!!! If there is an error, it is likely because the version of .msh file is above 2.
start_time = time.perf_counter()

msh = ReadGmsh(mesh_path)
msh = Mesh(msh)

print("ElementBoundary=", msh.GetBoundaries())

#%%
# Define the equation and BCs
fes = H1(msh, order=1)

u = fes.TrialFunction()
v = fes.TestFunction()

a = BilinearForm(fes, symmetric=True)
a += cuThermCond*grad(u)*grad(v)*dx
a += heatTransCoeff * u * v * ds("bc-pipe-htc")

f = LinearForm(fes)
f += surfHeatFlux * v * ds("bc-top-heatflux") + heatTransCoeff * T_inf * v * ds("bc-pipe-htc")

#%%
# Run the simulation

pre = preconditioners.Local(a)
a.Assemble()
f.Assemble()

gfu = GridFunction(fes)

ngsglobals.msg_level = 1
inv = CGSolver(a.mat, pre.mat, precision = 1e-10, printrates=True)
gfu.vec.data = inv*f.vec

end_time = time.perf_counter()
run_time = end_time - start_time
print(f'NGSolve run time = {run_time:.3f} seconds')

# #%%
# Save the results and compare with analytical solution

sol = gfu.vec.FV().NumPy()

res = pv.read(mesh_path)
points = res.points

res["sol"] = sol

res.save(output_path)