function C = ss_to_dss_state_space_transformed_default ...
    (eqn, opts, fac, opFac, B, opB)
%% function C = dss_to_ss_state_space_transformed_default...
%     (eqn, opts, fac, opFac, B, opB)
%
% This function returns either C = EL*B or C = EU*B, where EL and EU
% are the LU factors of E_ = EL*EU given in eqn.
% Matrices EL and EU are assumed to be quadratic.
%
% Inputs
%   eqn             struct contains data for equations
%
%   opts            struct contains parameters for the algorithm
%
%   fac             character specifying the used factor
%                       fac = 'L' performs op(EL)*op(B)
%                       fac = 'U' performs op(EU)*op(B)
%
%   opFac           character specifying the shape of the used factor
%                       opFac = 'N' performs EL*op(B) or EU*op(B)
%                       opFac = 'T' performs EL'*op(B) or EU'*op(B)
%
%   B               n-x-p matrix
%
%   opB             character specifying the shape of B
%                       opB = 'N' performs op(EL)*B or op(EU)*B
%                       opB = 'T' performs op(EL)*B' or op(EU)*B'
%
% Output
%   C = op(EL)*op(B) or op(EU)*op(B)
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
mess_assert(opts, ischar(fac) && ischar(opFac) && ischar(opB), ...
            'error_arguments', ...
            'fac, opFac or opB is not a char');

fac   = upper(fac);
opFac = upper(opFac);
opB   = upper(opB);

mess_assert(opts, (fac == 'L') || (fac == 'U'), ...
            'error_arguments', ...
            'fac is not ''N'' or ''T''');

mess_assert(opts, (opFac == 'N') || (opFac == 'T'), ...
            'error_arguments', ...
            'opFac is not ''N'' or ''T''');

mess_assert(opts, (opB == 'N') || (opB == 'T'), ...
            'error_arguments', ...
            'opB is not ''N'' or ''T''');

mess_assert(opts, isnumeric(B) && ismatrix(B), ...
            'error_arguments', ...
            'B has to ba a matrix');

%% Check data in eqn structure.
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

rowE = size_default(eqn, opts);
colE = rowE;

%% Perform solve operation.
if eqn.haveE % Case of non-identity E matrix.
    switch fac
        case 'L'
            switch opFac
                case 'N'
                    switch opB
                        case 'N' % Implement operation EL*B.
                            mess_assert(opts, colE == size(B, 1), ...
                                        'error_arguments', ...
                                        ['number of columns of E_ differs ' ...
                                         'with number of rows of B']);
                            C = eqn.ER(:, eqn.Ep) * (eqn.EL * B);

                        case 'T' % Implement operation EL*B'.
                            mess_assert(opts, colE == size(B, 2), ...
                                        'error_arguments', ...
                                        ['number of columns of E_ differs ' ...
                                         'with number of columns of B']);
                            C = eqn.ER(:, eqn.Ep) * (eqn.EL * B');
                    end

                case 'T'
                    switch opB
                        case 'N' % Implement operation EL'*B.
                            mess_assert(opts, rowE == size(B, 1), ...
                                        'error_arguments', ...
                                        ['number of rows of E_ differs with ' ...
                                         'number rows of B']);
                            C = eqn.EL' * (eqn.ER(:, eqn.Ep)' * B);
                        case 'T' % Implement operation EL'*B'.
                            mess_assert(opts, rowE == size(B, 2), ...
                                        'error_arguments', ...
                                        ['number of rows of E_ differs with ' ...
                                         'number of columns of B']);
                            C = eqn.EL' * (eqn.ER(:, eqn.Ep)' * B');
                    end

            end
        case 'U'
            switch opFac
                case 'N'
                    switch opB
                        case 'N' % Implement operation EU*B.
                            mess_assert(opts, colE == size(B, 1), ...
                                        'error_arguments', ...
                                        ['number of columns of E_ differs ' ...
                                         'with number of rows of B']);
                            C = eqn.EU * B(eqn.Eq, :);
                        case 'T' % Implement operation EU*B'.
                            mess_assert(opts, colE == size(B, 2), ...
                                        'error_arguments', ...
                                        ['number of columns of E_ differs ' ...
                                         'with number of columns of B']);
                            C = eqn.EU * B(:, eqn.Eq)';
                    end

                case 'T'
                    switch opB
                        case 'N' % Implement operation EU'*B.
                            mess_assert(opts, rowE == size(B, 1), ...
                                        'error_arguments', ...
                                        ['number of rows of E_ differs with ' ...
                                         'number rows of B']);
                            C(eqn.Eq, :) = eqn.EU' * B;
                        case 'T' % Implement operation EU'*B'.
                            mess_assert(opts, rowE == size(B, 2), ...
                                        'error_arguments', ...
                                        ['number of rows of E_ differs with ' ...
                                         'number of columns of B']);
                            C(eqn.Eq, :) = eqn.EU' * B';
                    end

            end
    end
else % Case of E_ = I_n, was set by init.
    C = B;
end
