M-M.E.S.S. - The Matrix Equation Sparse Solver Library for MATLAB and Octave
============================================================================

M-M.E.S.S. provides low rank solvers for large scale symmetric matrix
equations with sparse or sparse + low rank coefficients. The main
focus is on differential and algebraic Riccati equations appearing in
control and model order reduction, as well as algebraic Lyapunov
equations for, e.g., balanced truncation.

The underlying dynamical system may be of first or second order and
structured proper differential algebraic equations (DAEs) that allow
for implicit index reduction are also supported.

The solvers philosophy is to always work on the implicitly linearized
(for second order systems) and/or implicitly projected (in the DAE case)
matrix equations. That means the implicit Lyapunov or Riccati equation
is always of the form known for a standard first order ODE, that may
have a non identity but invertible E matrix.

See `help mess` for an overview of supported matrix equations and
system structures. 

Further, M-M.E.S.S. provides functions for Balanced Truncation and
(tangential) iterative rational Krylov algorithm (IRKA) for model order
reduction (MOR) of first order state space systems and some examples
demonstrate the use of the algorithms in MOR of second order systems and DAEs.

In close relation to the predecessor LyaPack, we use user supplied
functions that implement the actions of the system matrices E and A in
multiplication and (shifted) solves. We provide those functions for
standard state space systems, second order systems, structured DAEs of
index 1 and 2, as well as second order DAEs of index 1, 2 and 3.

Copyright 2009-2020 
 by Jens Saak, Martin Koehler, Peter Benner (MPI Magdeburg)

The software is licensed under GPLv2 or later. See [LICENSE.md](LICENSE.md)
and [COPYING](COPYING) for details.


## Installation Instructions

See [INSTALL.md](INSTALL.md) for details.

## Getting started

Change to the installation directory, run `mess_path` and check `help mess`
for the basic information about supported matrix equations and core solvers.

In case you need functionality beyond that of `mess_lyap` and `mess_care`,
consult the  demonstration routines in the DEMOS folder for example use
cases of the other and underlying solvers.

## Contact

[WWW](https://www.mpi-magdeburg.mpg.de/projects/mess), 
[GIT](https://gitlab.mpi-magdeburg.mpg.de/mess/mmess-releases)

email: [Jens Saak](mailto:saak@mpi-magdeburg.mpg.de), 
       [Martin Köhler](mailto:koehlerm@mpi-magdeburg.mpg.de)

## Citation

See [CITATION.md](CITATION.md) for details about citing the software.

## Further reading
- J. Saak, M. Voigt; '''Model reduction of constrained mechanical
  systems in M-M.E.S.S'''; ''IFAC-PapersOnLine'', 51(2):
  661--666, 2018. [https://doi.org/10.1016/j.ifacol.2018.03.112](https://doi.org/10.1016/j.ifacol.2018.03.112) 
