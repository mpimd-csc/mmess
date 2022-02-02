function X = sol_E_state_space_transformed_default(eqn, opts, opE, B, opB)
%% function X = sol_E_state_space_transformed_default(eqn,opts,opE,B,opB)
%
% This function returns X = E_\B, where matrix E_ given by structure eqn
% and input matrix B could be transposed.
% Matrix E_ is assumed to be quadratic.
%
% Inputs
%   eqn             struct contains data for equations
%
%   opts            struct contains parameters for the algorithm
%
%   opE             character specifying the shape of E_
%                       opE = 'N' solves E_*X = opB(B)
%                       opE = 'T' solves E_'*X = opB(B)
%
%   B               m-x-p matrix
%
%   opB             character specifying the shape of B
%                       opB = 'N' solves opE(E_)*X = B
%                       opB = 'T' solves opE(E_)*X = B'
%
% Output
%    X              matrix solving opE(E_)*X = opB(B)
%
% This function uses another default function size_default(eqn, opts) to
% obtain the number of rows of matrix E_ in structure eqn.

%
% This file is part of the M-M.E.S.S. project 
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%


%% Check input parameters.
assert(ischar(opE) && ischar(opB), ...
    'MESS:error_arguments', ...
    'opE or opB is not a char');

opE = upper(opE);
opB = upper(opB);

assert((opE == 'N') || (opE == 'T'), ...
    'MESS:error_arguments', ...
    'opE is not ''N'' or ''T''');

assert((opB == 'N') || (opB == 'T'), ...
    'MESS:error_arguments', ...
    'opB is not ''N'' or ''T''');

assert(isnumeric(B) && ismatrix(B), ...
    'MESS:error_arguments', ...
    'B has to ba a matrix');

%% Check data in eqn structure.
assert(isfield(eqn, 'EL'), ...
        'MESS:error_arguments', ...
        'field eqn.EL is not defined');
    assert(isfield(eqn, 'EU'), ...
        'MESS:error_arguments', ...
        'field eqn.EU is not defined');

rowE = size_default(eqn, opts);
colE = rowE;

%% Perform solve operation.
switch opE
    case 'N'
        switch opB
            case 'N' % Implement solve E_*X = B.
                assert(rowE == size(B, 1), ...
                    'MESS:error_arguments', ...
                    ['number of rows of E_ differs with ' ...
                    'number rows of B']);
                X = eqn.EU \ (eqn.EL \ B);

            case 'T' % Implement solve E_*X = B'.
                assert(rowE == size(B, 2), ...
                    'MESS:error_arguments', ...
                    ['number of rows of E_ differs with ' ...
                    'number of columns of B']);
                X = eqn.EU \ (eqn.EL \ B');
        end

    case 'T'
        switch opB
            case 'N' % Implement solve E_'*X = B.
                assert(colE == size(B, 1), ...
                    'MESS:error_arguments', ...
                    ['number of columns of E_ differs with ' ...
                    'number of rows of B']);
                X = eqn.EL' \ (eqn.EU' \ B);

            case 'T' % Implement solve E_'*X = B'.
                assert(colE == size(B, 2), ...
                    'MESS:error_arguments', ...
                    ['number of columns of E_ differs with ' ...
                    'number of columns of B']);
                X = eqn.EL' \ (eqn.EU' \ B');
        end

end
