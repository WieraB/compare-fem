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
case_name = "case2"
mesh_path = "./meshes/mesh_square.msh"

output_path = f"./output/{case_name}/{case_name}_ngsolve.vtu"

#%%
# Load the mesh
# WARNING!!! If there is an error, it is likely because the version of .msh file is above 2.
start_time = time.perf_counter()

msh = ReadGmsh(mesh_path)
msh = Mesh(msh)

print("ElementBoundary=", msh.GetBoundaries())

#%%
# Define the equation and BCs
fes = H1(msh, order=1, dirichlet="Left|Right")

u = fes.TrialFunction()
v = fes.TestFunction()

a = BilinearForm(fes, symmetric=True)
a += grad(u)*grad(v)*dx

f = LinearForm(fes)
g = sin(6*x)
f += 10 * g * v * ds("Top")
f += 10 * g * v * ds("Bottom")

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

#%%
# Save the results

sol = gfu.vec.FV().NumPy()

res = pv.read(mesh_path)
res["sol"] = sol

res.save(output_path)