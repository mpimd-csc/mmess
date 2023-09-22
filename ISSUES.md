# Known Issues

## Known problems with M-M.E.S.S. itself

* It is essential to remove all earlier M-M.E.S.S. versions from the
  search path and uninstall previous versions of the
  toolbox. Otherwise conflicts are almost guaranteed.
* The octave dense Riccati solvers can not handle indefinite
  right-hand-sides. These can occur in the
  `mess_solve_projected_eqn` routine.
* The `mess_care` has been observed to work less accurate in
  octave, when octave is compiled with GCC before version 5. That
  means, on RHEL/CentOS/ScientificLinux 6/7, Ubuntu 14.04, SLES 11/12
  one should not rely on the default system GCC for the octave
  compilation. This may also affect other routines.
* The `dae_2_so` and `dae_3_so` usfs compute the wrong residuals when
  `mess_res2_norms` is used. This can be fixed along the lines of the
  current `dae_2` implementation. Internal residual computations do
  not use `mess_res2_norms`, but rely on known residual factors.
* The current `dae_1` implementation is overwriting `eqn.B` and
  `eqn.C` with the ones on the hidden manifold.
* The `mess_splitting_dre` requires the input `eqn.Rinv`. This implements
 `E'*d/dt X(t)*E = -C'*C - E'*X(t)*A - A'*X(t)*E + E'*X(t)*B*Rinv*B'*X(t)*E (T)`.
  The input is not used for the computation of the feedback matrix
  `K(t) = E'*X(t)*B`.
* `mess_bdf_dre` and `mess_rosenbrock_dre.m` do not use the input
  `eqn.Rinv`.
* `mess_splitting_dre` does not give the expected order of
  approximation in the LTV case for order larger 2.
* `mess_splitting_dre` works for the LTV case only if
   `A(t)*A(s)=A(s)*A(t)`. This is not checked in the code and has to
   be ensured by the user.
* BDF solvers can currently not be used with line-search activated in
  the inner Newton-Kleinman solves.
* `mess_lrnm`in LDL_T mode with line-search is sometimes unstable in Octave.
* `mess_tangential_irka` was observed to converge exceptionally slow
  using MATLAB R2019b on certain Intel Sandy bridge processors.
* `mess_lrradi` can crash in rare cases inside
  `mess_RADI_get_shifts_hamOpti_generalized`, probably due to an
  economy size QR with incompatible dimensions.
* `mess_res2_norms` was observed to be less accurate
  using MATLAB R2019b on certain Intel Westmere processors.
* `exp_action` in the splitting methods does not work with its
  'Krylov' method when `dae_2` usfs are used.

## Compatibility with other Software

* The sssMOR toolbox from MORLab @ TUM contains an earlier version of
  M-M.E.S.S. this may lead to conflicts with our internal
  data-management. Also some wrong number of inputs or outputs error
  messages can appear.
