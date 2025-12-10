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
case_name = "case6"
mesh_path = "./meshes/mesh_square_QUAD8.msh"
output_path = f"./output/{case_name}/{case_name}_fenics.vtu"

#%%
# Define initial distribution and material properties
u_inf = 0.0
topDisp = 1.0e-3  # m


# Material properties
cuEMod= 108e9   # Pa
cuPRatio = 0.33     # -
mu = cuEMod / (2 * (1 + cuPRatio))
lambda_ = cuEMod * cuPRatio / ((1 - 2 * cuPRatio) * (1 + cuPRatio))

order = 2

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

V = fem.functionspace(msh, ("Lagrange", order, (msh.geometry.dim,)))

# left_tag = 1
# facets_left = facet_marker.find(left_tag)
# right_tag = 2
# facets_right = facet_marker.find(right_tag)
top_tag = 3
facets_top = facet_marker.find(top_tag)
bottom_tag = 4
facets_bottom = facet_marker.find(bottom_tag)

# dofs_left = fem.locate_dofs_topological(V=V, entity_dim=1, entities=facets_left)
# dofs_right = fem.locate_dofs_topological(V=V, entity_dim=1, entities=facets_right)
dofs_top = fem.locate_dofs_topological(V=V, entity_dim=1, entities=facets_top)
dofs_bottom = fem.locate_dofs_topological(V=V, entity_dim=1, entities=facets_bottom)


bc_bottom = fem.dirichletbc(value=np.array([0, 0], dtype=dolfinx.default_scalar_type), dofs=dofs_bottom, V=V)
bc_top = fem.dirichletbc(value=np.array([0, topDisp], dtype=dolfinx.default_scalar_type), dofs=dofs_top, V=V)

T = fem.Constant(msh, dolfinx.default_scalar_type((0, 0))) # Stress vector at the boundary, and can be prescribed as a BC

ds = ufl.Measure("ds", domain=msh)

# Strain function
def epsilon(u):
    return ufl.sym(
        ufl.grad(u)
    )  # Equivalent to 0.5*(ufl.nabla_grad(u) + ufl.nabla_grad(u).T)

# Stress function
def sigma(u):
    return lambda_ * ufl.nabla_div(u) * ufl.Identity(len(u)) + 2 * mu * epsilon(u)


v = ufl.TestFunction(V)
uh = fem.Function(V)
f = fem.Constant(msh, dolfinx.default_scalar_type((0, 0))) # Body force per unit volume (e.g. weight)
a = ufl.inner(sigma(uh), epsilon(v)) * ufl.dx
L = ufl.dot(f, v) * ufl.dx + ufl.dot(T, v) * ds

uh.x.array[:] = u_inf

F = a - L

#%%
# Run the simulation

petsc_options = {
    "snes_type": "newtonls",
    "snes_linesearch_type": "none",
    "snes_atol": 1e-6,
    "snes_rtol": 1e-6,
    "snes_monitor": None,
    "ksp_error_if_not_converged": True,
    "ksp_type": "cg",
    "ksp_rtol": 1e-6,
    "ksp_monitor": None,
    "pc_type": "hypre",
    "pc_hypre_type": "boomeramg",
    "pc_hypre_boomeramg_max_iter": 1,
    "pc_hypre_boomeramg_cycle_type": "v",
}

problem = NonlinearProblem(
    F,
    uh,
    bcs=[bc_bottom, bc_top],
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
# Calculate stresses

# W = fem.functionspace(msh, ("DG", 1, (msh.geometry.dim, msh.geometry.dim))) # Tensor-valued discontinuous space

W1 = fem.functionspace(msh, ("Lagrange", order, (msh.geometry.dim, msh.geometry.dim)))
W2 = fem.functionspace(msh, ("Lagrange", order))

sigma_h = fem.Function(W1)
stress_expr = sigma(uh)
projector = fem.Expression(stress_expr, W1.element.interpolation_points)
sigma_h.interpolate(projector)

s = sigma(uh) - 1.0 / 3 * ufl.tr(sigma(uh)) * ufl.Identity(len(uh))
von_Mises = ufl.sqrt(3.0 / 2 * ufl.inner(s, s))
stress_expr = fem.Expression(von_Mises, W2.element.interpolation_points)
stresses = fem.Function(W2)
stresses.interpolate(stress_expr)

#%%
# Save the results
# Extracting solution vector from solution function uh 
# and re-arranging dofs to match mesh dof order to match what we did in NGSolve example

# with io.XDMFFile(msh.comm, f"./output/{case_name}/{case_name}_fenics.xdmf", "w") as xdmf:
#     xdmf.write_mesh(msh)
#     uh.name = "Displacement"
#     xdmf.write_function(uh)
#     sigma_h.name = "Stress"
#     xdmf.write_function(sigma_h)

res = pv.read(mesh_path)
points = res.points

points_unordered = V.tabulate_dof_coordinates()
sol_unordered = uh.x.array.reshape((points_unordered.shape[0], 2))
print(sol_unordered.shape)

sol_unordered2 = sigma_h.x.array.reshape((-1, msh.geometry.dim, msh.geometry.dim))
print(sol_unordered2.shape)

sol_unordered3 = stresses.x.array
print(sol_unordered3.shape)

tree2 = KDTree(points_unordered)

_, indices = tree2.query(points)

sol = sol_unordered[indices]
sol2 = sol_unordered2[indices]
sol3 = sol_unordered3[indices]

res["disp"] = sol
res["stress"] = sol2
res["stress_von_Mises"] = sol3

res.save(output_path)