The **DEMOS/models/** directory contains the collection of benchmark
systems used in the demonstration examples.

**Data_Rail** The well known Rail aka Steel Profile example from the
  Oberwolfach Collection and a bilinear reformulation on the same
  geometry.

**FDM_2D** Contains the functions for generating scalable models of
  the finite difference semi-discretized model of a heat equation on
  the unit square. The functios support convection, as well as
  reaction terms at the users choice. They are exact copies of the
  files from the Lyapack package.

**NSE** Prepared to store the external download for the Karman vortex
  shedding model in a 2d channel.

**SingleChainMSD**
  A simple scalable mass spring damper system

**TripleCchain**
  The Truhar/Veselic model made from three coupled mass-spring-damper
  chains. Size, i.e. masses per chain, and damper viscosity, as well
  as parameters in the Rayleigh damping used here can be set by the
  user. An interesting parametrization can be found in
  **example_from_Saak09.m**

**msd_ind3_by_t_stykel**
  Tatjana Stykels mass spring damper system with holonomic
  constraints, i.e., the index-3 DAE case.

**stokes**
  The finite volume semidiscretized Stokes (index-2) model by Tatjana
  Stykel and Michael Schmidt.

The M.E.S.S. team would like to thank Tatjana Stykel for providing the
routines for the index-3 MSD model and the stokes model.
