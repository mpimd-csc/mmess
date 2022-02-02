# Dependencies of M-M.E.S.S. 2.2

## Basic requirements

MathWorks MATLAB R2014a and above, or GNU Octave 5.1 and above.

### Recommendation for MATLAB users

MATLAB R2017a has improved handling of negative definite matrices with in the
"backslash" operator. We recommend using this version or later ones for optimal
performance.

Note that MATLAB R2017a and R2017b also contain a bug in "backslash" that can
cause extraordinarily slow computations with certain block structured matrices.
Use

```
spparms('usema57', 0);
```

to fix this, or upgrade to at least R2017b update 5.

MATLAB R2021a has proven to be up to 20% faster than predecessors in
our continuous integration tests. So upgrading is highly recommended.

### Recommendations for Octave users

Some of the demonstration examples rely on MATLAB style ODE solvers
with proper support for non-trivial mass matrices. We highly recommend
Octave 5.2 and above with proper SUNDIALS support compiled in.

Octave 6.2.0 based on OpenBLAS has shown to be competitive with MATLAB
in recent releases for our continuous integration tests.

The splitting schemes further need Octave to be compiled with FFTW support.

## Optional dependencies

1.) Some functions can benefit from the Control Systems Toolbox in MATLAB or the
Control package in Octave. More precisely, whenever projected Lyapunov or
Riccati equations are solved, we use `lyap`, `care` or `icare` if available,
but fallbacks exist in case any of those is not found.
(affected files: helpers/mess_galerkin_projection_acceleration.m)

2.) Factorizations of symmetric and (semi)definite projected solutions in the
context of 1.) can benefit from the file *cholp.m* from the **Test
Matrix Computation Toolbox** by Nicholas J. Higham.
(affected files: helpers/mess_galerkin_projection_acceleration.m)
