function n = size_state_space_transformed_default(eqn, opts)
%% function n = size_default(eqn, opts)
%
% This function returns the number of rows of matrix A_ in structure eqn.
%
% Inputs
%   eqn             struct contains data for equations
%
%   opts            struct contains parameters for the algorithm
%
% Output
%    n              number of rows of matrix A_ in structure eqn
%
% This function does not use other default functions.

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

mess_assert(opts, isfield(eqn, 'A_'), ...
            'error_arguments', ...
            'field eqn.A_ is not defined');

n = size(eqn.A_, 1);
