# Known Issues 
* It is essential to remove all earlier M-M.E.S.S. versions from the
  search path and uninstall previous versions of the
  toolbox. Otherwise conflicts are almost guaranteed.
* The octave dense Riccati solvers can not handle indefinite
  right-hand-sides. These can occur in the
  `mess_galerkin_projection_acceleration` routine.
* The `mess_care` has been observed to work less accurate in
  octave, when octave is compile with GCC before version 5. That
  means, on RHEL/CentOS/ScientificLinux 6/7, Ubuntu 14.04, SLES 11/12
  one should not rely on the default system GCC for the octave
  compilation. This may also affect other routines.
* The `dae_2_so` and `dae_3_so` usfs compute the wrong residuals when
  `mess_res2_norms` is used. This can be fixed along the lines of the
  current `dae_2` implementation. Internal residual computations do
  not use `mess_res2_norms`, but rely on known residual factors.
* The current `dae_1` implementation is overwriting `eqn.B` and
  `eqn.C` with the ones on the hidden manifold.

# Compatibility with other Software

* The sssMOR toolbox from MORlab @ TUM contains an earlier version of
  M-M.E.S.S. this may lead to conflicts with our internal
  data-management. Also some wrong number of inputs or outputs error
  messages can appear.
  
  
