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
case_name = "case6"
mesh_path = "./meshes/mesh_square_QUAD8.msh"

output_path = f"./output/{case_name}/{case_name}_ngsolve.vtu"

#%%
# Define material properties, initial distribution, and solver parameters

u_inf = 0.0
topDisp = 1.0e-3  # m

# Material properties
cuEMod= 108e9   # Pa
cuPRatio = 0.33     # -
mu = cuEMod / (2 * (1 + cuPRatio))
lambda_ = cuEMod * cuPRatio / ((1 - 2 * cuPRatio) * (1 + cuPRatio))

# Solver parameters
maxits = 100
tol = 1e-2
precision = 1e-6
num_threads = 1

#%%
# Load the mesh
# WARNING!!! If there is an error, it is likely because the version of .msh file is above 2.
start_time = time.perf_counter()

msh = ReadGmsh(mesh_path)
msh = Mesh(msh)

print("ElementBoundary=", msh.GetBoundaries())

#%%
# Define the equation and BCs

def Stress(strain):
    return 2*mu*strain + lambda_*Trace(strain)*Id(2)    

fes = VectorH1(msh, order=2, dirichlet="bottom|top")
gfu = GridFunction(fes)
gfu.Set((u_inf, u_inf))
gfu.Set((0.0, 0.0), definedon=msh.Boundaries("bottom"))
gfu.Set((0.0, topDisp), definedon=msh.Boundaries("top"))
res = gfu.vec.CreateVector()
du = gfu.vec.CreateVector()

u = fes.TrialFunction()
v = fes.TestFunction()

#%%
# Run the simulation

SetNumThreads(num_threads)

with TaskManager():

    for it in range(maxits):
        print ("Iteration {:3}  ".format(it),end="")
    
        a = BilinearForm(InnerProduct(Stress(Sym(Grad(u))), Sym(Grad(v))).Compile()*dx)
        a.Assemble()
    
        # force = CF( (0.0, 0) )
        # f = LinearForm(force*v*ds("left"))
        # f.Assemble()
    
        pre = Preconditioner(a, type="local")
        inv = CGSolver(a.mat, pre.mat, precision = precision)
    
        # res.data = a.mat * gfu.vec - f.vec
        res.data = a.mat * gfu.vec
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
# n_dofs = gfu.dof_map()


print(sol.shape)

res = pv.read(mesh_path)
points = res.points

point = msh(0.,        1.6302521, 0.       )
point = msh(points[100, 0], points[100, 1], points[100, 2])
print(points[100, :])
ab = gfu(point)
print (ab)

sol = sol.reshape((-1, 2))
print(sol.shape)

res["sol"] = sol

res.save(output_path)