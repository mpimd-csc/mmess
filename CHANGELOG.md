## version 2.0
- New RADI iteration for AREs
- New splitting methods for autonomous DREs
- New splitting and BDF methods for non-autonomous DREs
- New operator manager only requires non-empty functions and replaces
  non-existent ones with a general `mess_do_nothing` function 
- improved Riccati iteration
- updated minimum required/recommended Matlab and octave versions
  (see `DEPENDENCIES.md`)
- unified function interfaces for top level calls
- unified handling of low rank updated operators. Now always A+UV' is
  used. (Note the sign of the update and the transposition in V) 
- major updates in the MOR routines
- some resturcturing in the opts structure. 
  * `opts.adi.shifts` has moved to `opts.shifts` such that also RADI
    can use it independent of ADI
  *  opts.norm now determines the norm for all methods rather than
     having to consistently specifiy the same norm in each substructure
  * initial feedbacks for the Riccati solvers are now stored in the
    `opts` structure for the method rather than `eqn` 
- The projection shift routine uses the flag `opts.shifts.implicitVtAV`. 
  Default is `true`. If set to `false` A*V is computed explicitly.
- several consistency updates and bug fixes
- general code cleaning and pretty printing
- redesign of the demos
  - turned scripts into actual demo functions
  - new demos for indefinite AREs and H-infinity control
- CI testing
  - demos serve as system tests 
  - additional unit tests for the smaller building blocks and backend routines

---
## version 1.0.1
- Minor consistency and bug fixes and improved integrity of metafiles.
- updated documentation 
- Removed replacements directory since its content was not needed for
  Matlab after release 2010b and Octave after 4.0. 

---
## version 1.0
Compared to the predecessor lyapack a couple of things have changed.
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
- A tangential IRKA implementation for none-DAE systems was added
