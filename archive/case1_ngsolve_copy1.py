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
case_name = "case1"
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
fes = H1(msh, order=1, dirichlet="Left|Right|Top|Bottom")

u = fes.TrialFunction()
v = fes.TestFunction()

a = BilinearForm(fes, symmetric=True)
a += grad(u)*grad(v)*dx

f = LinearForm(fes)
f += 32 * (y*(1-y)+x*(1-x)) * v * dx

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
# Save the results and compare with analytical solution

sol = gfu.vec.FV().NumPy()

res = pv.read(mesh_path)
points = res.points
x_coord = points[:, 0]
y_coord = points[:, 1]

sol_exact = 16*x_coord*(1-x_coord)*y_coord*(1-y_coord)
abs_error = np.abs(sol - sol_exact)
abs_error_max = np.max(abs_error)
abs_error_mean = np.mean(abs_error)
print(f"Max. absolute error between NGSolve and analytical solution : {abs_error_max:.3e}.")
print(f"Avg. absolute error between NGSolve and analytical solution : {abs_error_mean:.3e}.")

res["sol"] = sol
res['absolute_error'] = abs_error

res.save(output_path)