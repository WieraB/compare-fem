# compare-fem
Comparing the following FEM implementations - NGSolve, FEniCS, and MOOSE


## Cases


### Case 1 -- Poisson equation on a 2D domain, steady state, Dirichlet BCs and xy-dependent source, mesh_square.msh, serial, analytical solution available

* NGSolve
  * NGSolve run time = 0.024 seconds
  * Max. absolute error between NGSolve and analytical solution : 5.032e-04.
  * Avg. absolute error between NGSolve and analytical solution : 2.873e-05.
  * Max. absolute error between NGSolve and FEniCS : 9.271e-08.
  * Avg. absolute error between NGSolve and FEniCS : 5.592e-09.
  * Max. absolute error between NGSolve and MOOSE : 3.195e-07.
  * Avg. absolute error between NGSolve and MOOSE : 5.003e-08.

* FEniCS
  * FEniCS run time = 0.026 seconds
  * Max. absolute error between FEniCS and analytical solution : 5.031e-04.
  * Avg. absolute error between FEniCS and analytical solution : 2.872e-05.

* MOOSE
  * MOOSE run time = 2.298 seconds
  * Max. absolute error between MOOSE analytical solution : 5.031e-04.
  * Avg. absolute error between MOOSE and analytical solution : 2.873e-05.


### Case 2 -- Poisson equation on a 2D domain, steady state, Dirichlet and x dependent Neumann BCs, mesh_square.msh, serial

* NGSolve
  * NGSolve run time = 0.023 seconds
  * Max. absolute error between NGSolve and FEniCS : 6.022e-07.
  * Avg. absolute error between NGSolve and FEniCS : 1.207e-07.
  * Max. absolute error between NGSolve and MOOSE : 4.777e-11.
  * Avg. absolute error between NGSolve and MOOSE : 1.005e-11.

* FEniCS
  * FEniCS run time = 0.027 seconds

* MOOSE
  * MOOSE run time = 2.350 seconds

### Case 3 -- Heat conduction on a 3D domain, steady state, Neumann BCs, mesh_block_pipe_refined.msh, serial

* NGSolve
  * NGSolve run time = 0.255 seconds
  * Max. absolute error between NGSolve and FEniCS : 2.085e-05.
  * Avg. absolute error between NGSolve and FEniCS : 1.625e-06.
  * Max. absolute error between NGSolve and MOOSE : 9.006e-09.
  * Avg. absolute error between NGSolve and MOOSE : 2.216e-09.

* FEniCS
  * FEniCS run time = 1.190 seconds

* MOOSE
  * MOOSE run time = 3.389 seconds


