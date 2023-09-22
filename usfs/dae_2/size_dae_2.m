function n = size_dae_2(eqn, opts, oper) %#ok<INUSD>
% function n = size_dae_2(eqn, opts, oper)
%
% This function returns the number of rows of the implicitly projected A
% matrix of the index-2 system.
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
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

if not(isfield(eqn, 'manifold_dim'))    || not(isnumeric(eqn.manifold_dim))
    mess_err(opts, 'error_arguments', ...
             ['Missing or corrupted manifold_dim field detected in ' ...
              'equation structure.']);
end
n = eqn.manifold_dim;

end
