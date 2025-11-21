# solve the Poisson equation -Delta u = f
# with Dirichlet boundary condition u = 0

from ngsolve import *
# from netgen.geom2d import unit_square

ngsglobals.msg_level = 1

# generate a triangular mesh of mesh-size 0.2
mesh = Mesh(unit_square.GenerateMesh(maxh=0.2))

# H1-conforming finite element space
fes = H1(mesh, order=1, dirichlet=[1,2,3,4])

# define trial- and test-functions
u = fes.TrialFunction()
v = fes.TestFunction()

# the right hand side
f = LinearForm(fes)
f += 32 * (y*(1-y)+x*(1-x)) * v * dx

# the bilinear-form 
a = BilinearForm(fes, symmetric=True)
a += grad(u)*grad(v)*dx

a.Assemble()
f.Assemble()

# the solution field 
gfu = GridFunction(fes)
gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
# print (u.vec)


# plot the solution (netgen-gui only)
Draw (gfu)
Draw (-grad(gfu), mesh, "Flux")

exact = 16*x*(1-x)*y*(1-y)
print ("L2-error:", sqrt (Integrate ( (gfu-exact)*(gfu-exact), mesh)))

import numpy as np

sol = gfu.vec.FV().NumPy()

vertices = mesh.vertices
points = [vertex.point for vertex in vertices]
points = np.array(points)
x_coord = points[:, 0]
y_coord = points[:, 1]

sol_exact = 16*x_coord*(1-x_coord)*y_coord*(1-y_coord)
abs_error = np.abs(sol - sol_exact)
abs_error_max = np.max(abs_error)
abs_error_mean = np.mean(abs_error)
print(f"Max. absolute error between NGSolve and analytical solution : {abs_error_max:.3e}.")
print(f"Avg. absolute error between NGSolve and analytical solution : {abs_error_mean:.3e}.")

