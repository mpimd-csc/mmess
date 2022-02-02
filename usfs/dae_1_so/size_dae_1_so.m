function n = size_dae_1_so(eqn, opts, oper)%#ok<INUSD>
% function n = size_dae_1_so(eqn, opts, oper)
%
% This function returns the number of rows of matrix A in (2) in help
% mess_usfs_dae1_so.
%
%   Input:
%
%   eqn     struct contains data for equations
%
%   opts    struct contains parameters for the algorithm
%
%   oper    struct contains function handles for operation
%           with A and E
%
%   Output:
%
%   n       size of the implicitly projected A matrix
%
% This function does not use other dae_1_so functions.
%
% See also mess_usfs_dae_1_so


%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

if not(isfield(eqn, 'nd'))    || not(isnumeric(eqn.nd))
    error('MESS:nd',...
    'Missing or Corrupted nd field detected in equation structure.');
end
n = 2 * eqn.nd;
end
