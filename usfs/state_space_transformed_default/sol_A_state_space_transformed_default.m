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
%                   and if eqn.haveE == true
%                       opA = 'N' solves (EL\A_/EU)*X = opB(B)
%                       opA = 'T' solves (EU'\A_'/EL')*X = opB(B)
%
%   B               m-x-p matrix
%
%   opB             character specifying the shape of B
%                       opB = 'N' solves opA(A_)*X = B
%                       opB = 'T' solves opA(A_)*X = B'
%                   and if eqn.haveE == true
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
else
    eqn.haveE = false;
end

rowA = size_default(eqn, opts);
colA = rowA;

%% Perform solve operation.
if eqn.haveE
    switch opA
        case 'N'
            switch opB
                case 'N' % Implement solve (EL\A_/EU)*X = B.
                    mess_assert(opts, rowA == size(B, 1), ...
                                'error_arguments', ...
                                ['number of rows of A_ differs with ' ...
                                 'number rows of B']);

                    tempX(eqn.iEq(eqn.Aq), :) =  eqn.AU \ (eqn.AL \ ...
                                                           (eqn.AR(:, eqn.Ap) \ (eqn.ER(:, eqn.Ep) * ...
                                                                                 (eqn.EL * B))));
                    X = eqn.EU * tempX;

                case 'T' % Implement solve (EL\A_/EU)*X = B'.
                    mess_assert(opts, rowA == size(B, 2), ...
                                'error_arguments', ...
                                ['number of rows of A_ differs with ' ...
                                 'number of columns of B']);
                    tempX(eqn.iEq(eqn.Aq), :) =  eqn.AU \ (eqn.AL \ ...
                                                           (eqn.AR(:, eqn.Ap) \ (eqn.ER(:, eqn.Ep) * ...
                                                                                 (eqn.EL * B'))));
                    X = eqn.EU * tempX;
            end

        case 'T'
            switch opB
                case 'N' % Implement solve (EL\A_'/EU)'*X = B.
                    mess_assert(opts, colA == size(B, 1), ...
                                'error_arguments', ...
                                ['number of columns of A_ differs with ' ...
                                 'number of rows of B']);

                    tempX(eqn.iAq(eqn.Eq), :) = eqn.EU' * B;
                    X = eqn.EL' * (eqn.ER(:, eqn.Ep)' * (eqn.AR(:, eqn.Ap)' \ ...
                                                         (eqn.AL' \ (eqn.AU' \ tempX))));

                case 'T' % Implement solve (EL\A_'/EU)'*X = B'.
                    mess_assert(opts, colA == size(B, 2), ...
                                'error_arguments', ...
                                ['number of columns of A_ differs with ' ...
                                 'number of columns of B']);
                    tempX(eqn.iAq(eqn.Eq), :) = eqn.EU' * B';
                    X = eqn.EL' * (eqn.ER(:, eqn.Ep)' * (eqn.AR(:, eqn.Ap)' \ ...
                                                         (eqn.AL' \ (eqn.AU' \ tempX))));
            end

    end
else % No E so no transformation required
    switch opA
        case 'N'
            switch opB
                case 'N' % Implement solve A_*X = B.
                    mess_assert(opts, rowA == size(B, 1), ...
                                'error_arguments', ...
                                ['number of rows of A_ differs with ' ...
                                 'number rows of B']);

                    X(eqn.Aq, :) = eqn.AU \ (eqn.AL \ (eqn.AR(:, eqn.Ap) \ B));

                case 'T' % Implement solve A_*X = B'.
                    mess_assert(opts, rowA == size(B, 2), ...
                                'error_arguments', ...
                                ['number of rows of A_ differs with ' ...
                                 'number of columns of B']);

                    X(eqn.Aq, :) = eqn.AU \ (eqn.AL \ (eqn.AR(:, eqn.Ap) \ B'));

            end

        case 'T'
            switch opB
                case 'N' % Implement solve A_'*X = B.
                    mess_assert(opts, colA == size(B, 1), ...
                                'error_arguments', ...
                                ['number of columns of A_ differs with ' ...
                                 'number of rows of B']);

                    X = eqn.AR(:, eqn.Ap)' \ (eqn.AL' \ (eqn.AU' \ B(eqn.Aq, :)));

                case 'T' % Implement solve A_'*X = B'.
                    mess_assert(opts, colA == size(B, 2), ...
                                'error_arguments', ...
                                ['number of columns of A_ differs with ' ...
                                 'number of columns of B']);

                    X = eqn.AR(:, eqn.Ap)' \ (eqn.AL' \ (eqn.AU' \ B(:, eqn.Aq)'));

            end

    end
end
