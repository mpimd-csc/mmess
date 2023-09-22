# Contributors

## Core Team

### Scientific Advisor

- Peter Benner [ORCID:0000-0003-3362-4103](https://orcid.org/0000-0003-3362-4103)

### Core Developers

- Jens Saak [ORCID:0000-0001-5567-9637](https://orcid.org/0000-0001-5567-9637)
- Martin Köhler [ORCID:0000-0003-2338-9904](https://orcid.org/0000-0003-2338-9904)

## Version 3.0

- Quirin Aumann [ORCID:0000-0001-7942-5703](https://orcid.org/0000-0001-7942-5703)
   - extended IRKA tests.
   - logger fixes.
- Björn Baran [ORCID:0000-0001-6570-3653](https://orcid.org/0000-0001-6570-3653)
   - fixes and improvements in DRE methods,
- Christian Himpe [ORCID:0000-0003-2194-6754](https://orcid.org/0000-0003-2194-6754)
   - code review
   - documentation fixes
   - release testing
- Martin Köhler [ORCID:0000-0003-2338-9904](https://orcid.org/0000-0003-2338-9904)
   - code review
- Davide Palitta [ORCID:0000-0002-6987-4430](https://orcid.org/0000-0002-6987-4430)
   - Prototype KSM for LEs and AREs.
- Jens Saak [ORCID:0000-0001-5567-9637](https://orcid.org/0000-0001-5567-9637)
   - Revised sparse-dense Sylvester solvers and extended test routine.
   - Testing, revision and optimization of the KSM codes.
   - Testing, revision and optimization of the logger.
   - Testing and revision of the iterative usfs.
   - restructured CI setup.
   - Revision of MOR methods and analysis functions.
   - code style and quality improvements.
   - automated style checker.
   - documentation updates.
   - release management.
   - updated MOR functions
     (unified interface, stability updates, merged backends)
   - reduced code duplication
   - revised error codes
   - `mess_para` issue fixes
   - experimental `default_iter` and `so_iter` usfs; supervision and
     testing
   - revised demonstration examples and benchmark model fetching.
- Steffen Werner [ORCID:0000-0003-1667-4862](https://orcid.org/0000-0003-1667-4862)
   - LDL_T KSM.
   - code reviews.

### Student Assistants and Interns

- Sebastian Bresch
   - Sparse-dense Sylvester solvers and basic test routine.
- Ronald Mendez
   - basic support for iterative linear solver support, via
     `default_iter` and `so_iter`usfs.
- Adrian Schulze
   - new logger framework
   - automated spellchecking
   - extended code style CI testing
   - code style improvements

## Version 2.2

- Quirin Aumann[ORCID:0000-0001-7942-5703](https://orcid.org/0000-0001-7942-5703)
   - bug and documentation fixes in IRKA
   - release testing
- Christian Himpe [ORCID:0000-0003-2194-6754](https://orcid.org/0000-0003-2194-6754)
   - code review
   - documentation fixes
   - release testing
- Jens Saak [ORCID:0000-0001-5567-9637](https://orcid.org/0000-0001-5567-9637)
   - improved MOR functions,
   - larger Rail examples (for both the linear and bilinear cases),
   - BIPS example fixes,
   - Documentation updates,
   - release testing,
   - code review
- Tony Stillfjord [ORCID:0000-0001-6123-4271](https://orcid.org/0000-0001-6123-4271)
   - splitting scheme for DREs related improvements.
- Steffen Werner [ORCID:0000-0003-1667-4862](https://orcid.org/0000-0003-1667-4862)
   - code review

### Student Assistants and Interns

- Adrian Schulze
   - Spellchecker for comments, strings and MD-files, for both CLI and CI

## Version 2.1

- Björn Baran [ORCID:0000-0001-6570-3653](https://orcid.org/0000-0001-6570-3653)
   - fixes and improvements in DRE methods,
   - Newton and ADI
- Christian Bertram [ORCID:0000-0002-9227-4580](https://orcid.org/0000-0002-9227-4580)
   - performance improvement in RADI
- Christian Himpe [ORCID:0000-0003-2194-6754](https://orcid.org/0000-0003-2194-6754)
   - code review
   - documentation fixes
   - release testing
- Jens Saak [ORCID:0000-0001-5567-9637](https://orcid.org/0000-0001-5567-9637)
   - improved documentation,
   - improved user feedback,
   - bug fixes,
   - improved CI setup,
   - improved MOR functions,
   - rewritten `dae_1_so` usfs,
   - refactored rail demo model,
   - bilinear BT demo and performance optimization,
   - automated packaging.
- Tony Stillfjord [ORCID:0000-0001-6123-4271](https://orcid.org/0000-0001-6123-4271)
   - minor update in splitting schemes for DREs.
- Steffen Werner [ORCID:0000-0003-1667-4862](https://orcid.org/0000-0003-1667-4862)
   - LDL^T RADI,
   - fixed and new demo for unstable Riccati equations,
   - fixed initial solution bugs in RADI,
   - changed `mess_care` backend to RADI.

### Student Assistants and Interns

- Sebastian Bresch
   - low-rank "bilinear Lyapunov" aka "Lyapunov plus positive" equation
     solver,
   - basic sparss and mechss support,
   - new usfs CI test framework.
- Adrian Schulze
   - code coverage report generation,
   - improved runtime reporting,
   - automatic packaging system.

## Version 2.0.1

- Björn Baran [ORCID:0000-0001-6570-3653](https://orcid.org/0000-0001-6570-3653)
   - DRE method fixes.
- Christian Himpe [ORCID:0000-0003-2194-6754](https://orcid.org/0000-0003-2194-6754)
   - code review and documentation fixes.
- Jens Saak [ORCID:0000-0001-5567-9637](https://orcid.org/0000-0001-5567-9637)
   - improved MOR functions,
   - partial release automation.
- Steffen Werner [ORCID:0000-0003-1667-4862](https://orcid.org/0000-0003-1667-4862)
   - bug fix for DAE_1 usfs.

## Version 2.0

- Björn Baran [ORCID:0000-0001-6570-3653](https://orcid.org/0000-0001-6570-3653)
   - BDF methods for non-autonomous DREs,
   - system tests.
- Patrick Kuerschner [ORCID:0000-0002-6114-8821](https://orcid.org/0000-0002-6114-8821)
   - RADI.
- Jens Saak [ORCID:0000-0001-5567-9637](https://orcid.org/0000-0001-5567-9637)
   - improved MOR functions,
   - test framework,
   - unit and system tests,
   - code and toolbox restructuring.
- Tony Stillfjord [ORCID:0000-0001-6123-4271](https://orcid.org/0000-0001-6123-4271)
   - splitting schemes for DREs.
- Steffen Werner [ORCID:0000-0003-1667-4862](https://orcid.org/0000-0003-1667-4862)
   - RADI,
   - improved Operator Manager,
   - improved Riccati iteration.

## Version 1.0 & 1.0.1

### Student Assistants and Interns

- Björn Baran [ORCID:0000-0001-6570-3653](https://orcid.org/0000-0001-6570-3653)
   - LDL^T based Algorithms and Differential Equations.
- Maximilian Behr [ORCID:0000-0001-8519-1632](https://orcid.org/0000-0001-8519-1632)
   - Operator Manager,
   - DAE function handles.
- Manuela Hund [ORCID:0000-0003-2888-3717](https://orcid.org/0000-0003-2888-3717)
   - Documentation.
- Steffen Werner [ORCID:0000-0003-1667-4862](https://orcid.org/0000-0003-1667-4862)
   - Riccati Iteration.

### Indirect Contributions

- Patrick Kürschner [ORCID:0000-0002-6114-8821](https://orcid.org/0000-0002-6114-8821)
   - experimental prototype codes for:
      - adaptive shifts,
      - residual factor based algorithms,
      - non-symmetric equations,
      - RADI.
- Norman Lang  [ORCID:0000-0002-9074-0103](https://orcid.org/0000-0002-9074-0103)
   - experimental prototype codes for:
      - LDL^T based algorithms,
      - Differential Lyapunov and Riccati equations.
- Heiko Weichelt [ORCID:0000-0002-9074-0103](https://orcid.org/0000-0002-9074-0103)
   - experimental prototype codes for:
      - inexact Newton with line-search
