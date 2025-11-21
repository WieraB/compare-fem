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
case_name = "case1"
mesh_path = "./meshes/mesh_square.msh"
output_path = f"./output/{case_name}/{case_name}_fenics.vtu"

#%%
# Load the mesh
# WARNING!!! It's important to specify that the mesh is 2D, otherwise it will perceive it as 3D

start_time = time.perf_counter()

gmsh.initialize()
gmsh.clear()
gmsh.model.add("loaded_mesh")

gmsh.open(mesh_path)

msh_data = dolfinx.io.gmsh.model_to_mesh(gmsh.model, MPI.COMM_SELF, 0, gdim=2)
msh = msh_data.mesh
cell_marker = msh_data.cell_tags
facet_marker = msh_data.facet_tags

print(f"Unique facet markers: {np.unique(facet_marker.values)}")

gmsh.finalize()

#%%
# Define the equation and BCs

V = fem.functionspace(msh, ("Lagrange", 1))

left_tag = 1
facets_left = facet_marker.find(left_tag)
right_tag = 2
facets_right = facet_marker.find(right_tag)
top_tag = 3
facets_top = facet_marker.find(top_tag)
bottom_tag = 4
facets_bottom = facet_marker.find(bottom_tag)

dofs_left = fem.locate_dofs_topological(V=V, entity_dim=1, entities=facets_left)
dofs_right = fem.locate_dofs_topological(V=V, entity_dim=1, entities=facets_right)
dofs_top = fem.locate_dofs_topological(V=V, entity_dim=1, entities=facets_top)
dofs_bottom = fem.locate_dofs_topological(V=V, entity_dim=1, entities=facets_bottom)

bc_left = fem.dirichletbc(value=ScalarType(0), dofs=dofs_left, V=V)
bc_right = fem.dirichletbc(value=ScalarType(0), dofs=dofs_right, V=V)
bc_top = fem.dirichletbc(value=ScalarType(0), dofs=dofs_top, V=V)
bc_bottom = fem.dirichletbc(value=ScalarType(0), dofs=dofs_bottom, V=V)

v = ufl.TestFunction(V)

uh = fem.Function(V)

x = ufl.SpatialCoordinate(msh)
a = ufl.inner(ufl.grad(uh), ufl.grad(v)) * ufl.dx

f = 32 * (x[1]*(1-x[1])+x[0]*(1-x[0]))

L = ufl.inner(f, v) * ufl.dx

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
    bcs=[bc_left, bc_right, bc_top, bc_bottom],
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
# Save the results and compare with analytical solution
# Extracting solution vector from solution function uh 
# and re-arranging dofs to match mesh dof order to match what we did in NGSolve example

res = pv.read(mesh_path)
points = res.points
x_coord = points[:, 0]
y_coord = points[:, 1]

sol_unordered = uh.x.array
points_unordered = V.tabulate_dof_coordinates()

tree2 = KDTree(points_unordered)

_, indices = tree2.query(points)

sol = sol_unordered[indices]

sol_exact = 16*x_coord*(1-x_coord)*y_coord*(1-y_coord)
abs_error = np.abs(sol - sol_exact)
abs_error_max = np.max(abs_error)
abs_error_mean = np.mean(abs_error)
print(f"Max. absolute error between FEniCS and analytical solution : {abs_error_max:.3e}.")
print(f"Avg. absolute error between FEniCS and analytical solution : {abs_error_mean:.3e}.")

res["sol"] = sol
res['absolute_error'] = abs_error

res.save(output_path)