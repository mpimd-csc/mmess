# Contributing to the M-M.E.S.S. project

Please note that this project is released with a Contributor [Code of
Conduct](CODE_OF_CONDUCT.md). By participating in this project you
agree to abide by its terms. Note further, that we follow
[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-v2.0%20adopted-ff69b4.svg)](code_of_conduct.md)
and would like to emphasize the following additional examples of
unacceptable behavior:

* Publishing or applying others' algorithms, implementations, or
  research ideas present in the private development repository
  (e.g. in a WIP branch or merge request), without their explicit
  permission.

When contributing to this project, please note the following style
conventions that are enforced (partially strictly enforced by the CI
server).

Most importantly, reflect your changes in the
[Changelog](CHANGELOG.md) so we do not have to reconstruct it from the
GIT history.

## Attribution

When you have contributed code to M-M.E.S.S. you and the content of
your contribution will be mentioned in the project's [contributors
file](CONTRIBUTORS.md).  Contributions are grouped by release, so if
you have contributed code to multiple releases, you will be mentioned
for each of these releases.

If, for some reason, you do not wish to be mentioned in the
[contributors file](CONTRIBUTORS.md), please give us a short note.
Also note that we cannot give attributions to trivial changes, such as
fixing a typo in the documentation or correcting a very simple
bug. However, your changes including your authorship will always be
included in the git history of the project.

## M-M.E.S.S. Style Guide

In general we prefer descriptive identifiers that support easier
reading and understanding of the code. Other than that, follow the
guidelines below:

**Naming of functions:**

* to avoid shadowing of functions from other toolboxes/packages or
  the MATLAB core, all M-M.E.S.S. solver functions start with `mess_`
* demonstration routines can have names deviating from the above
  principle but should be descriptive enough such that a basic idea
  of their contents is evident.
* user supplied functions need to have a special naming scheme
  described in the documentation or found in the `operatormanager`
  function inside the `usfs` folder
* Routines that do not obey the above have to be put into a
  `private` folder.

**Code Style:**

* We try to follow the Mathworks recommendation for maximum line
  lengths where possible. In exceptional cases lines of at most 120
  characters are allowed. The latter is enforced by our CI settings,
  i.e. the tests will fail when a longer line is
  detected. Documentation strings should strictly keep the 80
  characters bound for readability reasons.
* any checked in code must pass the `checkcode` tests in MATLAB,
  i.e. you are not allowed to have red or orange marks in the MATLAB
  editor. However, when absolutely necessary those warnings may be
  suppressed. (E.g. in ADI the growth of the solution factor is
  avoidable only with significand implementation overhead)
* We always use `not()` rather than `~` for negation for easier
  reading. (Enforced by CI!)
* To support our abstraction and provide maximum flexibility to users,
  all functions that are not explicitly intended for matrices need
  to have `eqn`, `opts`, and `oper` as mandatory input and output
  arguments. For example `mess_lradi` is a core Lyapunov solver that
  is intended for maximum abstraction and needs these arguments,
  while `mess_lyap` is intended as a large-scale replacement for
  MATLAB's `lyap` and thus intended for matrices.
* For the same abstraction reason, in core solvers, all operations
  with the system matrices have to be implemented via the `usfs`
  system and access to `eqn.A_` or `eqn.E_` is prohibited. (see `help
  mess_usfs` for details on the `usfs` system and philosophy)
* The options structure `opts` should contain a separate
  substructure for each algorithm/solver-function and options are
  only allowed on the top level when they are absolutely necessary
  for consistent operation. For example `mess_lrnm` and `mess_lradi`
  have their options in `opts.nm` and `opts.adi`. On the other hand,
  `opts.LDL_T` or `opts.norm` used to decide global settings like
  the shape of the factorization computed, or the norm that should
  consistently be used for measuring progress are allowed on the top
  level.
* Avoid code duplication. If code blocks or variations thereof need to
  be repeated, consider to put them into a function and call it
  where necessary. If a helper function lacks functionality consider
  extending rather than doubling it.

## Maintainability

Our GIT repository provides maximum flexibility for forking, branching,
merging, and commenting on work done. For the sake of maintainability
follow these guidelines:

* Try to break down things to as small issues as possible.
* Work on a single issue per merge request or branch.
* branch and merge often, since long living branches tend to become
  messy to handle and synchronize with the master.
* Try to 'rebase' your work onto the master before mergeing or
  requesting a merge, if possible

## Test Framework

We use an extensive test system on our continuous integration (CI)
server to ensure that new features do not break old functionality and
all features work in both MATLAB and GNU Octave, both with and without
toolboxes/packages installed. Currently CI time is limited to 2 hours
and all tests should be finished in that time.

We distinguish between so called unit tests and system tests. Here,
unit tests are testing the smallest possible building blocks. As a
rule of thumb, a routine that does not call external functions, or
only private functions is a candidate for a unit test.

Larger routines that run a hierarchy of other stuff, such as our
demonstrator functions in the `DEMOS` folder should be used to perform
the much bigger system tests.

* All files in `DEMOS` should be accompanied by a system test
  calling them on an, ideally, small example. (See comment about the
  CI time above)
* All `usfs` sets should be checked by unit tests. Use non-symmetric
  systems and system matrices where possible!
* Consider to provide a unit test for every helper routine in your
  method.
