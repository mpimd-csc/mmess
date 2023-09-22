function X = sol_E_state_space_transformed_default(eqn, opts, opE, B, opB)
%% function X = sol_E_state_space_transformed_default(eqn,opts,opE,B,opB)
%
% This function returns X = B; the input matrix B could be transposed.
% The transformed matrix E is assumed to be the identity in this function
% set. A non-identity E_ may still be present in the eqn structure and will
% be used in the transformation.
%
% Inputs
%   eqn             struct contains data for equations
%
%   opts            struct contains parameters for the algorithm
%
%   opE             character specifying the shape of the
%                   transformed E.
%                   unused since the transformed E acts as an identity.
%                   (still needs to be provided for consistency)
%
%   B               m-x-p matrix
%
%   opB             character specifying the shape of B
%                       opB = 'N' sets X = B
%                       opB = 'T' sets X = B'
%
% Output
%    X              matrix solving X = opB(B)
%
% This function uses another state_space_transformed_default function;
% size_state_space_transformed_default(eqn, opts) to  obtain the number of
% rows of the transformed matrix E from structure eqn.

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%% Check input parameters.
mess_assert(opts, ischar(opE) && ischar(opB), ...
            'error_arguments', ...
            'opE or opB is not a char');

opE = upper(opE);
opB = upper(opB);

mess_assert(opts, (opE == 'N') || (opE == 'T'), ...
            'error_arguments', ...
            'opE is not ''N'' or ''T''');

mess_assert(opts, (opB == 'N') || (opB == 'T'), ...
            'error_arguments', ...
            'opB is not ''N'' or ''T''');

mess_assert(opts, isnumeric(B) && ismatrix(B), ...
            'error_arguments', ...
            'B has to ba a matrix');

rowE = size_default(eqn, opts);

%% Perform solve operation.

switch opB
    case 'N' % Implement solve E_*X = B.
        mess_assert(opts, rowE == size(B, 1), ...
                    'error_arguments', ...
                    ['number of rows of E_ differs with ' ...
                     'number rows of B']);
        X = B;

    case 'T' % Implement solve E_*X = B'.
        mess_assert(opts, rowE == size(B, 2), ...
                    'error_arguments', ...
                    ['number of rows of E_ differs with ' ...
                     'number of columns of B']);
        X = B';
end
