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
%                   and if eqn.haveE == true
%                       opA = 'N' performs EL\A_/EU*opB(B)
%                       opA = 'T' performs EU'\A_'/EL'*opB(B)
%
%   B               m-x-p matrix
%
%   opB             character specifying the shape of B
%                       opB = 'N' performs opA(A_)*B
%                       opB = 'T' performs opA(A_)*B'
%                   and if eqn.haveE == true
%                       opB = 'N' performs opA(EL\A_/EU)*B
%                       opB = 'T' performs opA(EL\A_/EU)*B'
%
% Output
%   C = opA(A_)*opB(B) or opA(EL\A_/EU)*opB(B)

% This function uses another default function size_default(eqn, opts) to
% obtain the number of rows of matrix A_ in structure eqn.

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%% Check input parameters.
mess_assert(opts, ischar(opA) && ischar(opB), ...
            'error_arguments', ...
            'opA or opB is not a char');

opA = upper(opA);
opB = upper(opB);

mess_assert(opts, (opA == 'N') || (opA == 'T'), ...
            'error_arguments', ...
            'opA is not ''N'' or ''T''');

mess_assert(opts, (opB == 'N') || (opB == 'T'), ...
            'error_arguments', ...
            'opB is not ''N'' or ''T''');

mess_assert(opts, isnumeric(B) && ismatrix(B), ...
            'error_arguments', ...
            'B has to ba a matrix');

%% Check data in eqn structure.
mess_assert(opts, isfield(eqn, 'A_'), ...
            'error_arguments', ...
            'field eqn.A_ is not defined');

if isfield(eqn, 'haveE') && eqn.haveE
    mess_assert(opts, isfield(eqn, 'EL'), ...
                'error_arguments', ...
                'field eqn.EL is not defined');
    mess_assert(opts, isfield(eqn, 'EU'), ...
                'error_arguments', ...
                'field eqn.EU is not defined');
    mess_assert(opts, isfield(eqn, 'Ep'), ...
                'error_arguments', ...
                'field eqn.Ep is not defined');
    mess_assert(opts, isfield(eqn, 'Eq'), ...
                'error_arguments', ...
                'field eqn.Eq is not defined');
    mess_assert(opts, isfield(eqn, 'ER'), ...
                'error_arguments', ...
                'field eqn.ER is not defined');
else
    eqn.haveE = false;
end

rowA = size_default(eqn, opts);
colA = rowA;

%% Perform multiplication.
if eqn.haveE
    switch opA
        case 'N'
            switch opB
                case 'N' % Implement operation (EL\A_/EU)*B.
                    mess_assert(opts, colA == size(B, 1), ...
                                'error_arguments', ...
                                ['number of columns of A_ differs with ' ...
                                 'number of rows of B']);
                    tempC(eqn.Eq, :) = eqn.EU \ B;
                    C = eqn.EL \ (eqn.ER(:, eqn.Ep) \ (eqn.A_ * tempC));

                case 'T' % Implement operation (EL\A_/EU)*B'.
                    mess_assert(opts, colA == size(B, 2), ...
                                'error_arguments', ...
                                ['number of columns of A_ differs with ' ...
                                 'number of columns of B']);
                    tempC(eqn.Eq, :) = eqn.EU \ B';
                    C = eqn.EL \ (eqn.ER(:, eqn.Ep) \ (eqn.A_ * tempC));
            end

        case 'T'
            switch opB
                case 'N' % Implement operation (EL\A_/EU)'*B.
                    mess_assert(opts, rowA == size(B, 1), ...
                                'error_arguments', ...
                                ['number of rows of A_ differs with ' ...
                                 'number rows of B']);

                    tempC = eqn.A_' * (eqn.ER(:, eqn.Ep)' \ (eqn.EL' \ B));
                    C = eqn.EU' \ tempC(eqn.Eq, :);

                case 'T' % Implement operation (EL\A_/EU)'*B'.
                    mess_assert(opts, rowA == size(B, 2), ...
                                'error_arguments', ...
                                ['number of rows of A_ differs with ' ...
                                 'number of columns of B']);
                    tempC = eqn.A_' * (eqn.ER(:, eqn.Ep)' \ (eqn.EL' \ B'));
                    C = eqn.EU' \ tempC(eqn.Eq, :);
            end

    end
else % No E and thus no transformation required
    switch opA
        case 'N'
            switch opB
                case 'N' % Implement operation A_*B.
                    mess_assert(opts, colA == size(B, 1), ...
                                'error_arguments', ...
                                ['number of columns of A_ differs with ' ...
                                 'number of rows of B']);
                    C = eqn.A_ * B;

                case 'T' % Implement operation A_*B'.
                    mess_assert(opts, colA == size(B, 2), ...
                                'error_arguments', ...
                                ['number of columns of A_ differs with ' ...
                                 'number of columns of B']);
                    C = eqn.A_ * B';
            end

        case 'T'
            switch opB
                case 'N' % Implement operation A_'*B.
                    mess_assert(opts, rowA == size(B, 1), ...
                                'error_arguments', ...
                                ['number of rows of A_ differs with ' ...
                                 'number rows of B']);
                    C = eqn.A_' * B;

                case 'T' % Implement operation A_'*B'.
                    mess_assert(opts, rowA == size(B, 2), ...
                                'error_arguments', ...
                                ['number of rows of A_ differs with ' ...
                                 'number of columns of B']);
                    C = eqn.A_' * B';
            end

    end
end
