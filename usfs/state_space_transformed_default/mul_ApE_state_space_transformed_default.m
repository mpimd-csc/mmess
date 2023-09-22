function C = mul_ApE_state_space_transformed_default ...
    (eqn, opts, opA, p, opE, B, opB)
%% function C = mul_ApE_state_space_transformed_default ...
%    (eqn, opts, opA, p, opE, B, opB)
%
% This function returns C = EL\(A_ + pE_)/EU*B, where matrices A_ and E_
% given by structure eqn and input matrix B could be transposed.
% Matrices A_ and E_ are assumed to be quadratic.
%
% Inputs
%   eqn             struct contains data for equations
%
%   opts            struct contains parameters for the algorithm
%
%   opA             character specifying the shape of A_
%                   opA = 'N' performs EL\(A_ + p*op(E_))/EU * opB(B)
%                   opA = 'T' performs EL\(A_' + p*op(E_))/EU * opB(B)
%
%   p               scalar value
%
%   opE             character specifying the shape of E_
%                   opE = 'N' solves EL\(op(A_) + p*E_)/EU * opB(B)
%                   opE = 'T' solves EU'\(op(A_) + p*E_')/EL' * opB(B)
%
%   B               m-x-p matrix
%
%   opB             character specifying the shape of B
%                   opB = 'N' performs EL\(op(A_) + p*op(E_))/EU * B
%                   opB = 'T' performs EL\(op(A_) + p*op(E_))/EU * B'
%
% Output
%   C = EL\(op(A_) + p*op(E_))/EU * op(B)
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
mess_assert(opts, ischar(opA) && ischar(opE) && ischar(opB), ...
            'error_arguments', ...
            'opA or opB is not a char');

opA = upper(opA);
opE = upper(opE);
opB = upper(opB);

mess_assert(opts, (opA == 'N') || (opA == 'T'), ...
            'error_arguments', ...
            'opA is not ''N'' or ''T''');

mess_assert(opts, (opE == 'N') || (opE == 'T'), ...
            'error_arguments', ...
            'opE is not ''N'' or ''T''');

mess_assert(opts, (opB == 'N') || (opB == 'T'), ...
            'error_arguments', ...
            'opB is not ''N'' or ''T''');

mess_assert(opts, isnumeric(p) && (length(p) == 1), ...
            'error_arguments', ...
            'p is not a numeric scalar');

mess_assert(opts, isnumeric(B) && ismatrix(B), ...
            'error_arguments', ...
            'B has to ba a matrix');

%% Check data in eqn structure.
mess_assert(opts, isfield(eqn, 'A_'), ...
            'error_arguments', ...
            'field eqn.A_ is not defined');

mess_assert(opts, isfield(eqn, 'E_'), ...
            'error_arguments', ...
            'field eqn.E_ is not defined');

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

%% Perform multiplication.
if eqn.haveE % Case of non-identity E matrix.
    switch opA
        case 'N'
            switch opE
                case 'N'
                    switch opB
                        case 'N' % Implement EL\(A_ + pE_)/EU*B.
                            mess_assert(opts, colA == size(B, 1), ...
                                        'error_arguments', ...
                                        ['number of columns of A_ differs ' ...
                                         'with number of rows of B']);
                            tempC(eqn.Eq, :) = eqn.EU \ B;
                            C = eqn.EL \ (eqn.ER(:, eqn.Ep) \ ...
                                          ((eqn.A_ + p * eqn.E_) * tempC));
                        case 'T' % Implement EL\(A_ + pE_)/EU*B'.
                            mess_assert(opts, colA == size(B, 2), ...
                                        'error_arguments', ...
                                        ['number of columns of A_ differs ' ...
                                         'with number of columns of B']);
                            tempC(eqn.Eq, :) = eqn.EU \ B';
                            C = eqn.EL \ (eqn.ER(:, eqn.Ep) \ ...
                                          ((eqn.A_ + p * eqn.E_) * tempC));
                    end

                case 'T'

                    mess_err(opts, 'missing_feature', ...
                             ['The cases where opA differs from opE ' ...
                              'have not yet been implemented']);

            end
        case 'T'
            switch opE
                case 'N'

                    mess_err(opts, 'missing_feature', ...
                             ['The cases where opA differs from opE ' ...
                              'have not yet been implemented']);

                case 'T'
                    switch opB
                        case 'N' % Implement EU'\(A_' + pE_')/EL'*B.
                            mess_assert(opts, rowA == size(B, 1), ...
                                        'error_arguments', ...
                                        ['number of rows of A_ differs with ' ...
                                         'number rows of B']);
                            tempC = (eqn.A_' + p * eqn.E_') * (eqn.ER(:, eqn.Ep)' \ ...
                                                               (eqn.EL' \ B));
                            C = eqn.EU' \ tempC(eqn.Eq, :);

                        case 'T' % Implement EU'\(A_' + pE_')/EL'*B'.
                            mess_assert(opts, rowA == size(B, 2), ...
                                        'error_arguments', ...
                                        ['number of rows of A_ differs with ' ...
                                         'number of columns of B']);
                            tempC = (eqn.A_' + p * eqn.E_') * (eqn.ER(:, eqn.Ep)' \ ...
                                                               (eqn.EL' \ B'));
                            C = eqn.EU' \ tempC(eqn.Eq, :);
                    end

            end
    end
else % Case of E_ = I_n uses eqn.I_ set by init and does not need
    % sate space transformation
    switch opA
        case 'N'
            switch opB
                case 'N' % Implement (A_ + pI_n)*B.
                    mess_assert(opts, colA == size(B, 1), ...
                                'error_arguments', ...
                                ['number of columns of A_ differs with ' ...
                                 'number of rows of B']);
                    C = (eqn.A_ + p * eqn.I_) * B;
                case 'T' % Implement (A_ + pI_n)*B'.
                    mess_assert(opts, colA == size(B, 2), ...
                                'error_arguments', ...
                                ['number of columns of A_ differs with ' ...
                                 'number of columns of B']);
                    C = (eqn.A_ + p * eqn.I_) * B';
            end

        case 'T'
            switch opB
                case 'N' % Implement (A_' + pE_)*B.
                    mess_assert(opts, rowA == size(B, 1), ...
                                'error_arguments', ...
                                ['number of rows of A_ differs with ' ...
                                 'number rows of B']);
                    C = (eqn.A_' + p * eqn.I_) * B;
                case 'T' % Implement (A_' + pE_)*B'.
                    mess_assert(opts, rowA == size(B, 2), ...
                                'error_arguments', ...
                                ['number of rows of A_ differs with ' ...
                                 'number of columns of B']);
                    C = (eqn.A_' + p * eqn.I_) * B';
            end

    end
end

end
