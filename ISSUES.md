# Known Issues 

* The octave dense Riccati solvers can not handle indefinite
  right-hand-sides. These can occur in the
  `mess_galerkin_projection_acceleration` routine.
* The `dae_2_so` and `dae_3_so` usfs compute the wrong residuals when
  `mess_res2_norms` is used. This can be fixed along the lines of the
  current `dae_2` implementation.
* The current `dae_1` implementation is overwriting `eqn.B` and
  `eqn.C` with the ones on the hidden manifold.
