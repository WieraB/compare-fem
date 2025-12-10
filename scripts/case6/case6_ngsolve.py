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
cuEMod = 108e9   # Pa
cuPRatio = 0.33     # -
mu = cuEMod / (2 * (1 + cuPRatio))
lambda_ = cuEMod * cuPRatio / ((1 - 2 * cuPRatio) * (1 + cuPRatio))

# Solver parameters
maxits = 100
tol = 1e-2
precision = 1e-6
num_threads = 1
dim = 2

#%%
# Load the mesh
# WARNING!!! If there is an error, it is likely because the version of .msh file is above 2.
start_time = time.perf_counter()

msh = ReadGmsh(mesh_path)
msh = Mesh(msh)

print("ElementBoundary=", msh.GetBoundaries())

#%%
# Define the equation and BCs

# Strain function
def epsilon(u):
    return Sym(Grad(u))

# Stress function
def sigma(u):
    return lambda_ * Trace(epsilon(u)) * Id(dim) + 2 * mu * epsilon(u) 

fes = VectorH1(msh, order=dim, dirichlet="bottom|top")
print(fes.ndof/dim)
gfu = GridFunction(fes)
gfu.Set((u_inf, u_inf))
gfu.Set((0.0, 0.0), definedon=msh.Boundaries("bottom"))
gfu.Set((0.0, topDisp), definedon=msh.Boundaries("top"))


# gfu.components[0].Set(0.0, definedon=msh.Boundaries("bottom"))
# gfu.components[1].Set(0.0, definedon=msh.Boundaries("bottom"))
# gfu.components[0].Set(0.0, definedon=msh.Boundaries("top"))
# gfu.components[1].Set(topDisp, definedon=msh.Boundaries("top"))

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
    
        a = BilinearForm(fes, symmetric=True)
        a += InnerProduct(sigma(u), epsilon(v)).Compile()*dx
        a.Assemble()
    
        force = CF( (0.0, 0.0) )
        T = CF( (0.0, 0.0) )
        f = LinearForm(fes)
        f += force * v * dx + T * v * ds
        f.Assemble()
    
        pre = Preconditioner(a, type="local")
        inv = CGSolver(a.mat, pre.mat, precision = precision)
        # solvers.BVP(bf=a, lf=f, gf=gfu, pre=pre)
    
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

#%%
# Save the results

res = pv.read(mesh_path)
points = res.points

sol = np.zeros((points.shape[0], 2))

for i in range(points.shape[0]):
    point = msh(points[i, 0], points[i, 1], points[i, 2])
    sol[i, :] = gfu(point)

res["disp"] = sol

res.save(output_path)
