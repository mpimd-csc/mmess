function n = size_so_iter(eqn, opts, oper) %#ok<INUSD>
% function n = size_so_iter(eqn, opts, oper)
%
% This function returns the number of rows of matrix A_ in structure eqn.
%
%   Input:
%
%   eqn     structure contains data for equations
%
%   opts    structure contains parameters for the algorithm
%
%   oper    structure contains function handles for operation
%           with A and E
%
%   Output:
%
%   n      number of rows of matrix A_ in structure eqn
%
% This function does not use other so_iter functions.

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

n = size(eqn.M_, 1);
end
