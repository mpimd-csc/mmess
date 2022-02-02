% Function handles for solving state-space transformed matrix equations
% via Krylov subspace methods.
% Assume the linear first-order system of the form
%
%    Ex' = Ax + Bu,
%     y  = Cx,
%
% with E = LU invertible.
% Then in the state-space transformed case, we are considering the system
%
%    z' = (L\A/U)z + (L\B)u,
%    y  = (C/U)z.
%
% The system is encoded in the eqn structure
%
% The fieldnames for A and E have to end with _  to indicate that the data
% are inputdata for the algorithm.
%
% eqn.A_ = A
% eqn.E_ = E
% eqn.B  = B
% eqn.C  = C
%
% Note that E_ and A_ are expected to be sparse and of size n x n,
% while B and C may be dense (and will be converted to dense by
% some routines anyway) and should have far less columns (for B)
% and rows (for C) than n.
%

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%
