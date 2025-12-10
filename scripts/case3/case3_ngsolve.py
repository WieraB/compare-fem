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

# Solver parameters
maxits = 100
tol = 1e-6
precision = 1e-6

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
    a += cuThermCond*grad(u)*grad(v)*dx
    a += heatTransCoeff * u * v * ds("bc-pipe-htc")
    a.Assemble()

    f = LinearForm(fes)
    f += surfHeatFlux * v * ds("bc-top-heatflux") + heatTransCoeff * T_inf * v * ds("bc-pipe-htc")
    f.Assemble()

    pre = Preconditioner(a, type="local")

    ngsglobals.msg_level = 1
    inv = CGSolver(a.mat, pre.mat, precision = precision,printrates=True)

    res.data = a.mat * gfu.vec - f.vec
    res_norm = np.linalg.norm(res.FV().NumPy())
    if it == 0:
            res0_norm = res_norm

    du.data = inv * (-res)
    gfu.vec.data += du
    
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

# #%%
# Save the results

sol = gfu.vec.FV().NumPy()

res = pv.read(mesh_path)
points = res.points

res["sol"] = sol

res.save(output_path)