function C = mul_E_state_space_transformed_default(eqn, opts, opE, B, opB)
%% function C = mul_E_state_space_transformed_default(eqn,opts,opE,B,opB)
%
% This function returns C = B, input matrix B could be transposed.
% Matrix E_ may exist in the eqn structure, but is ignored since the
% transformed system is assumed to be in standard state space form, i.e.
% the transformed E is the identity.
%
% Inputs
%   eqn             struct contains data for equations
%
%   opts            struct contains parameters for the algorithm
%
%   opE             character specifying the transposition of the
%                   transformed E.
%                   unused since the transformed E acts as an identity.
%                   (still needs to be provided for consistency)
%
%   B               m-x-p matrix
%
%   opB             character specifying the shape of B
%                       opB = 'N' sets C=B
%                       opB = 'T' sets C=B'
%
% Output
%   C = opB(B)
%
% This function uses another user supplied function
% (size_state_space_transformed_default(eqn, opts) to obtain the number of
% rows of the transformed E matrix.

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

%% Check data in eqn structure.

rowE = size_state_space_transformed_default(eqn, opts);
colE = rowE;

%% Perform multiplication.
switch opB
    case 'N' % Implement operation E_*B.
        mess_assert(opts, colE == size(B, 1), ...
                    'error_arguments', ...
                    ['number of columns of E_ differs with ' ...
                     'number of rows of B']);
        C = B;

    case 'T' % Implement operation E_*B'.
        mess_assert(opts, colE == size(B, 2), ...
                    'error_arguments', ...
                    ['number of columns of E_ differs with ' ...
                     'number of columns of B']);
        C = B';

end
