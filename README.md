# compare-fem
Comparing the following FEM implementations - NGSolve, FEniCS, and MOOSE


## Cases


### Case 1 -- Linear Poisson equation on a 2D domain, steady state, Dirichlet BCs and xy-dependent source, mesh_square.msh, serial, analytical solution available

* NGSolve
  * NGSolve run time = 0.025 seconds
  * Max. absolute error between NGSolve and FEniCS : 1.505e-07.
  * Avg. absolute error between NGSolve and FEniCS : 2.236e-08.
  * Max. absolute error between NGSolve and MOOSE : 1.204e-07.
  * Avg. absolute error between NGSolve and MOOSE : 1.311e-08.
  * Max. absolute error between NGSolve and analytical solution : 5.031e-04.
  * Avg. absolute error between NGSolve and analytical solution : 2.872e-05.

* FEniCS
  * FEniCS run time = 0.023 seconds

* MOOSE
  * MOOSE run time = 2.279 seconds


### Case 2 -- Linear Poisson equation on a 2D domain, steady state, Dirichlet and x dependent Neumann BCs, mesh_square.msh, serial

* NGSolve
  * NGSolve run time = 0.029 seconds
  * Max. absolute error between NGSolve and FEniCS : 7.020e-07.
  * Avg. absolute error between NGSolve and FEniCS : 1.215e-07.
  * Max. absolute error between NGSolve and MOOSE : 1.750e-07.
  * Avg. absolute error between NGSolve and MOOSE : 3.641e-08.

* FEniCS
  * FEniCS run time = 0.025 seconds

* MOOSE
  * MOOSE run time = 2.281 seconds


### Case 3 -- Linear heat conduction on a 3D domain, steady state, Neumann BCs, mesh_block_pipe_refined.msh, serial (based on case09 in pyvale simcases)

* NGSolve
  * NGSolve run time = 0.483 seconds
  * Max. absolute error between NGSolve and FEniCS : 2.085e-05.
  * Avg. absolute error between NGSolve and FEniCS : 1.625e-06.
  * Max. absolute error between NGSolve and MOOSE : 5.476e-05.
  * Avg. absolute error between NGSolve and MOOSE : 1.221e-05.

* FEniCS
  * FEniCS run time = 1.013 seconds

* MOOSE
  * MOOSE run time = 3.291 seconds

### Case 4 -- Nonlinear heat conduction on a 3D domain, steady state, Neumann BCs, mesh_block_pipe_refined.msh, serial -- Nonlienar version of Case 3

* NGSolve
  * NGSolve run time = 2.582 seconds
  * Max. absolute error between NGSolve and FEniCS : 1.123e-04.
  * Avg. absolute error between NGSolve and FEniCS : 6.455e-05.
  * Max. absolute error between NGSolve and MOOSE : 2.382e-04.
  * Avg. absolute error between NGSolve and MOOSE : 2.211e-04.

* FEniCS
  * FEniCS run time = 2.122 seconds

* MOOSE
  * MOOSE run time = 10.690 seconds

### Case 5 -- Nonlinear heat conduction on a 3D domain, steady state, Neumann BCs, mesh_block_pipe_refined.msh, parallel (4 threads) -- Parallel version of Case 4

* NGSolve
  * NGSolve run time = 1.014 seconds
  * Max. absolute error between NGSolve and FEniCS : 1.122e-04.
  * Avg. absolute error between NGSolve and FEniCS : 6.455e-05.
  * Max. absolute error between NGSolve and MOOSE : 2.382e-04.
  * Avg. absolute error between NGSolve and MOOSE : 2.211e-04.

* FEniCS
  * FEniCS run time = 1.095 seconds

* MOOSE
  * MOOSE run time = 9.544 seconds

