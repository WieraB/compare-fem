from mpi4py import MPI
from petsc4py.PETSc import ScalarType
import numpy as np
import ufl
from dolfinx import fem, io, mesh, plot
from dolfinx.fem.petsc import LinearProblem, NonlinearProblem
import dolfinx
import meshio
import gmsh
import time
import numpy as np
from scipy.spatial import KDTree
import pyvista as pv

#%%
# Set file locations
case_name = "case3"
mesh_path = "./meshes/mesh_block_pipe_refined.msh"
output_path = f"./output/{case_name}/{case_name}_fenics.vtu"

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
# WARNING!!! It's important to specify that the mesh is 2D, otherwise it will perceive it as 3D

start_time = time.perf_counter()

gmsh.initialize()
gmsh.clear()
gmsh.model.add("loaded_mesh")

gmsh.open(mesh_path)

msh_data = dolfinx.io.gmsh.model_to_mesh(gmsh.model, MPI.COMM_SELF, 0, gdim=3)
msh = msh_data.mesh
cell_marker = msh_data.cell_tags
facet_marker = msh_data.facet_tags

print(f"Unique facet markers: {np.unique(facet_marker.values)}")

gmsh.finalize()

#%%
# Define the equation and BCs

V = fem.functionspace(msh, ("Lagrange", 1))

bc_top_heatflux_tag = 3
bc_pipe_htc_tag = 4

v = ufl.TestFunction(V)

uh = fem.Function(V)

ds = ufl.Measure("ds", domain=msh, subdomain_data=facet_marker)
a = cuThermCond * ufl.inner(ufl.grad(uh), ufl.grad(v)) * ufl.dx
a += heatTransCoeff * ufl.inner(uh, v) * ds(bc_pipe_htc_tag)

L = surfHeatFlux * v * ds(bc_top_heatflux_tag) + heatTransCoeff * T_inf * v * ds(bc_pipe_htc_tag)

F = a - L

#%%
# Run the simulation

petsc_options = {
    "snes_type": "newtonls",
    "snes_linesearch_type": "none",
    "snes_atol": 1e-10,
    "snes_rtol": 1e-10,
    "snes_monitor": None,
    "ksp_error_if_not_converged": True,
    "ksp_type": "cg",
    "ksp_rtol": 1e-10,
    "ksp_monitor": None,
    "pc_type": "hypre",
    "pc_hypre_type": "boomeramg",
    "pc_hypre_boomeramg_max_iter": 1,
    "pc_hypre_boomeramg_cycle_type": "v",
}

problem = NonlinearProblem(
    F,
    uh,
    bcs=[],
    petsc_options=petsc_options,
    petsc_options_prefix="nonlinpoisson",
)

problem.solve()
converged = problem.solver.getConvergedReason()
num_iter = problem.solver.getIterationNumber()
assert converged > 0, f"Solver did not converge, got {converged}."
print(
    f"Solver converged after {num_iter} iterations with converged reason {converged}."
)

end_time = time.perf_counter()
run_time = end_time - start_time
print(f'FEniCS run time = {run_time:.3f} seconds')

#%%
# Save the results
# Extracting solution vector from solution function uh 
# and re-arranging dofs to match mesh dof order to match what we did in NGSolve example

res = pv.read(mesh_path)
points = res.points

sol_unordered = uh.x.array
points_unordered = V.tabulate_dof_coordinates()

tree2 = KDTree(points_unordered)

_, indices = tree2.query(points)

sol = sol_unordered[indices]

res["sol"] = sol

res.save(output_path)