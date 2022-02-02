function C = mul_A_state_space_transformed_default(eqn, opts, opA, B, opB)
%% function C = mul_A_state_space_transformed_default(eqn,opts,opA,B,opB)
%
% This function returns C = A_*B, where matrix A_ given by structure eqn
% and input matrix B could be transposed.
% Matrix A_ is assumed to be quadratic.
%
% Inputs
%   eqn             struct contains data for equations
%
%   opts            struct contains parameters for the algorithm
%
%   opA             character specifying the shape of A_
%                       opA = 'N' performs A_*opB(B)
%                       opA = 'T' performs A_'*opB(B)
%                   and if eqn.haveE == 1
%                       opA = 'N' performs EL\A_/EU*opB(B)
%                       opA = 'T' performs EU'\A_'/EL'*opB(B)
%
%   B               m-x-p matrix
%
%   opB             character specifying the shape of B
%                       opB = 'N' performs opA(A_)*B
%                       opB = 'T' performs opA(A_)*B'
%                   and if eqn.haveE == 1
%                       opB = 'N' performs opA(EL\A_/EU)*B
%                       opB = 'T' performs opA(EL\A_/EU)*B'
%
% Output
%   C = opA(A_)*opB(B) or opA(EL\A_/EU)*opB(B)
%
% This function uses another default function size_default(eqn, opts) to
% obtain the number of rows of matrix A_ in structure eqn.

%
% This file is part of the M-M.E.S.S. project 
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
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

%% Perform multiplication.
if eqn.haveE
    switch opA
        case 'N'
            switch opB
                case 'N' % Implement operation (EL\A_/EU)*B.
                    assert(colA == size(B, 1), ...
                        'MESS:error_arguments', ...
                        ['number of columns of A_ differs with ' ...
                        'number of rows of B']);
                    C = eqn.EL \ (eqn.A_ * (eqn.EU \ B));

                case 'T' % Implement operation (EL\A_/EU)*B'.
                    assert(colA == size(B, 2), ...
                        'MESS:error_arguments', ...
                        ['number of columns of A_ differs with ' ...
                        'number of columns of B']);
                    C = eqn.EL \ (eqn.A_ * (eqn.EU \ B'));
            end

        case 'T'
            switch opB
                case 'N' % Implement operation (EL\A_/EU)'*B.
                    assert(rowA == size(B, 1), ...
                        'MESS:error_arguments', ...
                        ['number of rows of A_ differs with ' ...
                        'number rows of B']);
                    C = eqn.EU' \ (eqn.A_' * (eqn.EL' \ B));

                case 'T' % Implement operation (EL\A_/EU)'*B'.
                    assert(rowA == size(B, 2), ...
                        'MESS:error_arguments', ...
                        ['number of rows of A_ differs with ' ...
                        'number of columns of B']);
                    C = eqn.EU' \ (eqn.A_' * (eqn.EL' \ B'));
            end

    end
else
    switch opA
        case 'N'
            switch opB
                case 'N' % Implement operation A_*B.
                    assert(colA == size(B, 1), ...
                        'MESS:error_arguments', ...
                        ['number of columns of A_ differs with ' ...
                        'number of rows of B']);
                    C = eqn.A_ * B;

                case 'T' % Implement operation A_*B'.
                    assert(colA == size(B, 2), ...
                        'MESS:error_arguments', ...
                        ['number of columns of A_ differs with ' ...
                        'number of columns of B']);
                    C = eqn.A_ * B';
            end

        case 'T'
            switch opB
                case 'N' % Implement operation A_'*B.
                    assert(rowA == size(B, 1), ...
                        'MESS:error_arguments', ...
                        ['number of rows of A_ differs with ' ...
                        'number rows of B']);
                    C = eqn.A_' * B;

                case 'T' % Implement operation A_'*B'.
                    assert(rowA == size(B, 2), ...
                        'MESS:error_arguments', ...
                        ['number of rows of A_ differs with ' ...
                        'number of columns of B']);
                    C = eqn.A_' * B';
            end

    end
end
