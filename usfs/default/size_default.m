function n = size_default(eqn, opts, oper)%#ok<INUSD>
% function n = size_default(eqn, opts, oper)
%
% This function returns the number of rows of matrix A_ in structure eqn.
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
%   n       number of rows of matrix A_ in structure eqn
%
% This function does not use other default functions.

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2021 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

if(not(isfield(eqn,'A_')))
    error('MESS:error_arguments','field eqn.A_ is not defined');
end

n = size(eqn.A_,1);
end

