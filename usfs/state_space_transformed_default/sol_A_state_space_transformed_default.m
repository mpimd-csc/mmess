function X = sol_A_state_space_transformed_default(eqn, opts, opA, B, opB)
%% function X = sol_A_state_space_transformed_default(eqn,opts,opA,B,opB)
%
% This function returns X = A_\B, where matrix A_ given by structure eqn
% and input matrix B could be transposed.
% Matrix A_ is assumed to be quadratic.
%
% Inputs
%   eqn             struct contains data for equations
%
%   opts            struct contains parameters for the algorithm
%
%   opA             character specifying the shape of A_
%                       opA = 'N' solves A_*X = opB(B)
%                       opA = 'T' solves A_'*X = opB(B)
%                   and if eqn.haveE == 1
%                       opA = 'N' solves (EL\A_/EU)*X = opB(B)
%                       opA = 'T' solves (EU'\A_'/EL')*X = opB(B)
%
%   B               m-x-p matrix
%
%   opB             character specifying the shape of B
%                       opB = 'N' solves opA(A_)*X = B
%                       opB = 'T' solves opA(A_)*X = B'
%                   and if eqn.haveE == 1
%                       opB = 'N' solves opA(EL\A_/EU)*X = B
%                       opB = 'T' solves opA(EL\A_/EU)*X = B'
%
%    Output:
%
%    X              matrix solving opA(A_)*X = opB(B)
%                   or opA(EL\A_/EU)*X = op(B)
%
% This function uses another default function size_default(eqn, opts) to
% obtain the number of rows of matrix A_ in structure eqn.

%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, see <http://www.gnu.org/licenses/>.
%
% Copyright (C) Jens Saak, Martin Koehler, Peter Benner and others
%               2009-2020
%

%% Check input parameters.
assert(ischar(opA) && ischar(opB), ...
    'MESS:error_arguments', ...
    'opA or opB is not a char');

opA = upper(opA);
opB = upper(opB);

assert((opA == 'N') || (opA == 'T'), ...
    'MESS:error_arguments', ...
    'opA is not ''N'' or ''T''');

assert((opB == 'N') || (opB == 'T'), ...
    'MESS:error_arguments', ...
    'opB is not ''N'' or ''T''');

assert(isnumeric(B) && ismatrix(B), ...
    'MESS:error_arguments', ...
    'B has to ba a matrix');

%% Check data in eqn structure.
assert(isfield(eqn,'A_'), ...
    'MESS:error_arguments', ...
    'field eqn.A_ is not defined');

if isfield(eqn, 'haveE') && eqn.haveE
    assert(isfield(eqn, 'EL'), ...
        'MESS:error_arguments', ...
        'field eqn.EL is not defined');
    assert(isfield(eqn, 'EU'), ...
        'MESS:error_arguments', ...
        'field eqn.EU is not defined');
else
    eqn.haveE = 0;
end

rowA = size_default(eqn, opts);
colA = rowA;

%% Perform solve operation.
if eqn.haveE
    switch opA
        case 'N'
            switch opB
                case 'N' % Implement solve (EL\A_/EU)*X = B.
                    assert(rowA == size(B, 1), ...
                        'MESS:error_arguments', ...
                        ['number of rows of A_ differs with ' ...
                        'number rows of B']);
                    X = eqn.EU * (eqn.A_ \ (eqn.EL * B));
                    
                case 'T' % Implement solve (EL\A_/EU)*X = B'.
                    assert(rowA == size(B, 2), ...
                        'MESS:error_arguments', ...
                        ['number of rows of A_ differs with ' ...
                        'number of columns of B']);
                    X = eqn.EU * (eqn.A_ \ (eqn.EL * B'));
            end
            
        case 'T'
            switch opB
                case 'N' % Implement solve (EL\A_'/EU)'*X = B.
                    assert(colA == size(B, 1), ...
                        'MESS:error_arguments', ...
                        ['number of columns of A_ differs with ' ...
                        'number of rows of B']);
                    X = eqn.EL' * (eqn.A_' \ (eqn.EU' * B));
                    
                case 'T' % Implement solve (EL\A_'/EU)'*X = B'.
                    assert(colA == size(B, 2), ...
                        'MESS:error_arguments', ...
                        ['number of columns of A_ differs with ' ...
                        'number of columns of B']);
                    X = eqn.EL' * (eqn.A_' \ (eqn.EU' * B'));
            end
            
    end
else
    switch opA
        case 'N'
            switch opB
                case 'N' % Implement solve A_*X = B.
                    assert(rowA == size(B, 1), ...
                        'MESS:error_arguments', ...
                        ['number of rows of A_ differs with ' ...
                        'number rows of B']);
                    X = eqn.A_ \ B;
                    
                case 'T' % Implement solve A_*X = B'.
                    assert(rowA == size(B, 2), ...
                        'MESS:error_arguments', ...
                        ['number of rows of A_ differs with ' ...
                        'number of columns of B']);
                    X = eqn.A_ \ B';
            end
            
        case 'T'
            switch opB
                case 'N' % Implement solve A_'*X = B.
                    assert(colA == size(B, 1), ...
                        'MESS:error_arguments', ...
                        ['number of columns of A_ differs with ' ...
                        'number of rows of B']);
                    X = eqn.A_' \ B;
                    
                case 'T' % Implement solve A_'*X = B'.
                    assert(colA == size(B, 2), ...
                        'MESS:error_arguments', ...
                        ['number of columns of A_ differs with ' ...
                        'number of columns of B']);
                    X = eqn.A_' \ B';
            end
            
    end
end
