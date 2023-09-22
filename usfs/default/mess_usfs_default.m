% The first order system
%
%      E * z'(t) = A * z(t) + B * u(t)
%           y(t) = C * z(t)
%
% is encoded in the eqn structure
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
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%
