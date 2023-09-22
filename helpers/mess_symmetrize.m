function [A] = mess_symmetrize(A)
% MESS_SYMMETRIZE  makes sure the matrix A is numerically symmetric
%
% Input / Output
%
%  A  a theoretically symmetric matrix that may be numerically unsymmetric
%     and is symmetrized, i.e. numerically symmetric on output.

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

A = 0.5 * (A + A');
end
