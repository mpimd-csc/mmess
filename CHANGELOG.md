# CHANGELOG

## version 2.2

### Added

- Rail data are now available in additional sizes, via automatic
  downloads at runtime.

### Changed

- `mess_balanced_truncation` now adapts the default Lyapunov residual
  tolerances to the right hand side (RHS) and densifies RHS factors
  prior to compression.
- Documentation has been updated especially for usfs sets

### Fixed

- `mess_balanced_truncation` would fail on certain incomplete `opts`
  structures, which are now reporting proper errors so users can
  adapt.
- various smaller issues (outdated API, typographic errors) in the
  documentation have been fixed.
- automatic downloading of additional examples now also works in
  Octave
- `mess_get_BIPS` was failing on some models that used a differing
  naming scheme
- some comments in the `default` usfs have been updated to better match the
  current implementation
- figures now use consistent line-width for all graphs
- typographic errors and outdated defaults in the help texts have been
  corrected

## version 2.1

### Added

- Contribution guide & Code of Conduct
- README file for the usfs available via `help mess_usfs`.
- New solver for "Lyapunov plus positive" equations and related
  balancing based MOR of bilinear systems.
- basic support for sparss and mechss system classes.

### Changed

- We changed the license to BSD 2-Clause!
- `README.md` received some improvements pointing to references and
  documentation more explicitly.
- Sigma magnitude plots now allow presampling of the full order model
  when, e.g. a sequence of ROMs is being processed.
- `mess_balanced_truncation` now always performs a column compression
  of the right hand side factors in the two Gramian equations before
  starting the ADI.
- `mess_lradi` does not require precomputed shifts anymore. It
  computes them automatically independent of the type now.
- The RADI now allows for LDL^T factored right-hand sides in the Riccati
  equation and solution factors are expanded more efficiently.
- `LQR_FDM_unstable` is now `LQG_FDM_unstable_nwt` and got fixed for
  appropriate computation of dual regulator and filter Riccati equations in
  ZZ^T and LDL^T formats. In this context, the new demo
  `LQG_FDM_unstable_radi` has been added to show the same computations
  using the RADI method.
- Backend of `mess_care` changed from Newton method to RADI.
- time recording in the demos is safer now.
- dae_1_so usfs got rewritten completely and are documented much
  better now.
- documentation of the second order usfs was improved.
- the shift change criterion's default bound in tangential IRKA got
  lowered to 1e-6 (was 1e-2).

### Fixed

- `mess_balanced_truncation` would fail on certain incomplete `opts`
  structures, which are now filled automatically.
- The TV2 model in the second order index-3 demo file had a typo.
- `LICENSE.md` was showing a different license (fixed GPLv2) than the
  other files (GPLV2+)
- dead code was removed or received better documentation of its
  purpose.
- DRE solver received various minor fixes.
- the tangential IRKA received some minor fixes, such that IRKA
  options are used more reliably and shift change monitory is more correct.
- bugs in `mess_lrradi` for initial solutions and residuals in case of
  filter Riccati equation ('N') have been fixed.

## version 2.0.1

### Changed

- many function headers and help texts got improved/completed

### Fixed

- DAE_1 usfs failed for certain systems with non-symmetric A.
- LTV BDF could break in certain situations and was not following the
  general naming scheme for some variables.
- mess_res2_norms would break when more than 4 output arguments were requested

## version 2.0

### Added

- New RADI iteration for AREs
- New splitting methods for autonomous DREs
- New splitting and BDF methods for non-autonomous DREs
- New operator manager only requires non-empty functions and replaces
  non-existent ones with a general `mess_do_nothing` function
- renamed `opts.bdf.stage` to `opts.bdf.step`.
- CI testing
  - demos serve as system tests
  - additional unit tests for the smaller building blocks and backend routines

### Changed

- improved Riccati iteration
- updated minimum required/recommended Matlab and Octave versions
  (see `DEPENDENCIES.md`)
- unified function interfaces for top level calls
- unified handling of low rank updated operators. Now always A+UV' is
  used. (Note the sign of the update and the transposition in V)
- major updates in the MOR routines
- some restructuring in the opts structure.
  - `opts.adi.shifts` has moved to `opts.shifts` such that also RADI
    can use it independent of ADI
  - `opts.norm` now determines the norm for all methods rather than
     having to consistently specify the same norm in each substructure
  - initial feedbacks for the Riccati solvers are now stored in the
    `opts` structure for the method rather than `eqn`
- The projection shift routine uses the flag `opts.shifts.implicitVtAV`.
  Default is `true`. If set to `false` A*V is computed explicitly.
- redesign of the demos
  - turned scripts into actual demo functions
  - new demos for indefinite AREs and H-infinity control

### Fixed

- several consistency updates and bug fixes
- general code cleaning and pretty printing

## version 1.0.1

### Changed

- updated documentation
- Removed replacements directory since its content was not needed for
  Matlab after release 2010b and Octave after 4.0.

### Fixed

- Minor consistency and bug fixes and improved integrity of metafiles.
- CI testing
  - demos serve as system tests
  - additional unit tests for the smaller building blocks and backend routines

## version 1.0

Compared to the predecessor LyaPack a couple of things have changed.

- The user supplied functions are now managed by an operator manager
- The low rank ADI now has:
  - optimized treatment of E matrices in generalized equations
  - more choices for shift selection, including completely automatic
    generation of shifts
  - improved stopping criteria based on low rank factors of the current residual
  - automatic generation of real low rank factors also for complex shifts
- The Newton-Kleinman iteration features:
  - optimized treatment of E matrices in generalized equations
  - improved stopping criteria based on low rank factors of the current residual
  - inexact Newton, line search and Galerkin projection acceleration
- Examples have been extended
- The Riccati iteration for H-infinity Riccati equations was added
- DSPMR has not yet been ported to the new infrastructure
- The SRM routine for balanced truncation is only available for
  none-DAE systems. Still, DAE versions are included in the
  corresponding DEMOS.
- A tangential IRKA implementation for non-DAE systems was added

