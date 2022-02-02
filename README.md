# M-M.E.S.S. - The Matrix Equation Sparse Solver Library for MATLAB and Octave

M-M.E.S.S. provides low-rank solvers for large-scale symmetric matrix
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
functions (usfs) that implement the actions of the system matrices E and A in
multiplication and (shifted) solves. We provide those functions for
standard state space systems, second order systems, structured DAEs of
index 1 and 2, as well as second order DAEs of index 1, 2 and 3. For
more information on usfs see `help mess_usfs`.

Copyright 2009-2022
 by Jens Saak, Martin Koehler, Peter Benner (MPI Magdeburg)

The software uses a BSD 2-Clause license. See [LICENSE.md](LICENSE.md)
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
[GITLAB](https://gitlab.mpi-magdeburg.mpg.de/mess/mmess-releases),
[GITHUB](https://github.com/mpimd-csc/mmess)

[email](mailto:mess@mpi-magdeburg.mpg.de)

## Citation

See [CITATION.md](CITATION.md) for details about citing the software.

## Further reading

- P. Benner, M. Koehler, J. Saak, **Matrix equations, sparse solvers:
  M-M.E.S.S.-2.0.1 – philosophy, features and application for (parametric)
  model order reduction**,
  in: P. Benner, T. Breiten, H. Faßbender, M. Hinze, T. Stykel, R. Zimmermann
  (Eds.), *Model Reduction of Complex Dynamical Systems*, Vol. 171 of
  International Series of Numerical Mathematics, Birkhäuser, Cham, 2021,
  pp. 369–392.
  [https://doi.org/10.1007/978-3-030-72983-7_18](https://doi.org/10.1007/978-3-030-72983-7_18).
- J. Saak, M. Voigt, **Model reduction of constrained mechanical systems in
  M-M.E.S.S.**, *IFAC-PapersOnLine 9th Vienna International Conference on
  Mathematical Modelling MATHMOD 2018*,
  Vienna, Austria, 21–23 February 2018 51 (2) (2018) 661–666.
  [https://doi.org/10.1016/j.ifacol.2018.03.112](https://doi.org/10.1016/j.ifacol.2018.03.112).
- P. Benner, J. Saak, **Efficient solution of large scale Lyapunov and Riccati
  equations arising in model order reduction problems**,
  *Proc. Appl. Math. Mech.* 8 (1) (2008) 10085–10088.
  [https://doi.org/10.1002/pamm.200810085](https://doi.org/10.1002/pamm.200810085).

## History

- **2000 LyaPack**: M-M.E.S.S. originates in the work of Penzl and
  especially his software package LyaPack.
- **2003-2007** LyaPack 1.1 - 1.8 authored by Jens Saak improve the
  handling of non-identity E matrices.
- **2008** the first conference talk about the new project labeled
  M.E.S.S. is held at GAMM 2008 in Bremen (Germany).
- **2016 M-M.E.S.S.-1.0 and 1.0.1** first public releases of the
  greatly rewritten toolbox.
- **2019 M-M.E.S.S.-2.0** adds differential Riccati equations.
- **2020 M-M.E.S.S.-2.0.1** fixes several bugs and adds
  improvements for MOR.
- **2021 M-M.E.S.S.-2.1** adds Lyapunov plus positive equations and
  BT of bilinear systems.
- **2022 M-M.E.S.S.-2.2** fixes several smaller bugs and adds
  improvements to code style and performance, and improves documentation

## Roadmap

- **M-M.E.S.S.-3.0**
  - bilinear control problems
  - Krylov-projection-based solvers
  - sparse-dense Sylvester equations
- **M-M.E.S.S.-4.0**
  - sparse Sylvester equations
  - non-symmetric AREs

