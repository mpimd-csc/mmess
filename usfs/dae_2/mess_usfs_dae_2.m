% The first order system
%
%      E * z'(t) = A * z(t) + B * u(t)
%           y(t) = C * z(t)
%
% is encoded in the eqn structure
%
% The fieldnames for A and E have to end with _  to indicate that the data
% are inputdata for the algorithm. Further A_ and E_ have to be
% substructured as given below.
%
% eqn.A_ = [ A11 A12;
%            A21  0 ]
% eqn.E_ = [ E1  0;
%             0  0 ]
% eqn.B  = B
% eqn.C  = C
%
% The sizes of A11 and E1 have to coincide and the value needs to
% be specified in eqn.manifold_dim. Also B has eqn.manifold_dim rows and C eqn.manifold_dim
% columns.
% Furthermore, A12 needs to have full column-rank and A21 full row-rank.
%
% See the DEMOS/DAE2 folder for usage examples.
%

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%
