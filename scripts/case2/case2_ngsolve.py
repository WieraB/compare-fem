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
# Define initial distribution and solver parameters
T_inf = 0.0

# Solver parameters
maxits = 100
tol = 6.0e-01
precision = 1e-6

#%%
# Define the equation and BCs
fes = H1(msh, order=1, dirichlet="Left|Right")

gfu = GridFunction(fes)
gfu.Set(T_inf)
res = gfu.vec.CreateVector()
du = gfu.vec.CreateVector()

u = fes.TrialFunction()
v = fes.TestFunction()


#%%
# Run the simulation

for it in range(maxits):
    print ("Iteration {:3}  ".format(it),end="")

    a = BilinearForm(fes, symmetric=True)
    a += grad(u)*grad(v)*dx
    a.Assemble()

    f = LinearForm(fes)
    g = sin(6*x)
    f += 10 * g * v * ds("Top")
    f += 10 * g * v * ds("Bottom")
    f.Assemble()

    pre = Preconditioner(a, type="local")

    ngsglobals.msg_level = 1
    inv = CGSolver(a.mat, pre.mat, precision = precision,printrates=True)

    res.data = a.mat * gfu.vec - f.vec
    res_norm = np.linalg.norm(res.FV().NumPy())
    if it == 0:
            res0_norm = res_norm

    du.data = inv * (-res)
    gfu.vec.data += du  # Newton-like update
    
    rel_res = res_norm / res0_norm
    print(f"Relative residual = {rel_res:.2e}")

    if rel_res < tol:
        print("Converged after {} iterations".format(it + 1))
        break

if it == maxits - 1:
    print("Did not converge within {} iterations".format(maxits))
    print("Final relative residual = {:.2e}".format(rel_res))


end_time = time.perf_counter()
run_time = end_time - start_time
print(f'NGSolve run time = {run_time:.3f} seconds')

#%%
# Save the results

sol = gfu.vec.FV().NumPy()

res = pv.read(mesh_path)
res["sol"] = sol

res.save(output_path)