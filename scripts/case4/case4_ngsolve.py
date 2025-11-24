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
from scipy.interpolate import interp1d

#%%
# Set file locations
case_name = "case4"
mesh_path = "./meshes/mesh_block_pipe_refined.msh"

output_path = f"./output/{case_name}/{case_name}_ngsolve.vtu"

#%%
# Define material properties, loads, and solver parameters

# Thermal Loads/BCs
coolantTemp = 100.0      # degC
# heatTransCoeff = lambda temperature: np.exp(temperature/100) + 100
heatTransCoeff = lambda temperature: temperature * 100 + 100
surfHeatFlux = 5.0e6    # W.m^-2
T_inf = coolantTemp

# Material Properties: Pure (OFHC) Copper at 250degC
cuDensity = 8829.0  # kg.m^-3
cuThermCond = 384.0 # W.m^-1.K^-1
cuSpecHeat = 406.0  # J.kg^-1.K^-1

# Solver parameters
maxits = 100
tol = 1e-6

#%%
# Define heat transfer coefficient function

temps = np.arange(0, 3000, 0.5)
h_values = heatTransCoeff(temps)
h_func = interp1d(temps, h_values, fill_value="extrapolate")

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
h_gf = GridFunction(H1(msh,order=1))
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
    temp = gfu.vec.data
    h = h_func(temp)
    h_gf.vec.data = h
    h_cf = CoefficientFunction(h_gf)

    a = BilinearForm(fes, symmetric=True)
    a += cuThermCond*grad(u)*grad(v)*dx
    a += h_cf * u * v * ds("bc-pipe-htc")
    a.Assemble()

    f = LinearForm(fes)
    f += surfHeatFlux * v * ds("bc-top-heatflux") + h_cf * T_inf * v * ds("bc-pipe-htc")
    f.Assemble()

    pre = Preconditioner(a, type="local")
    inv = CGSolver(a.mat, pre.mat)

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

# #%%
# Save the results and compare with analytical solution

sol = gfu.vec.FV().NumPy()

res = pv.read(mesh_path)
points = res.points

res["sol"] = sol

res.save(output_path)