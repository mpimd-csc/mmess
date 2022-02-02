function C = ss_to_dss_state_space_transformed_default...
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
%                   	opFac = 'T' performs EL'*op(B) or EU'*op(B)
%
%   B               n-x-p matrix
%
%   opB             character specifying the shape of B
%                       opB = 'N' performs op(EL)*B or op(EU)*B
%                       opB = 'T' performs op(EL)*B' or op(EU)*B'
%
% Output
%	C = op(EL)*op(B) or op(EU)*op(B)
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
assert(ischar(fac) && ischar(opFac) && ischar(opB), ...
    'MESS:error_arguments', ...
    'fac, opFac or opB is not a char');

fac   = upper(fac);
opFac = upper(opFac);
opB   = upper(opB);

assert((fac == 'L') || (fac == 'U'), ...
    'MESS:error_arguments', ...
    'fac is not ''N'' or ''T''');

assert((opFac == 'N') || (opFac == 'T'), ...
    'MESS:error_arguments', ...
    'opFac is not ''N'' or ''T''');

assert((opB == 'N') || (opB == 'T'), ...
    'MESS:error_arguments', ...
    'opB is not ''N'' or ''T''');

assert(isnumeric(B) && ismatrix(B), ...
    'MESS:error_arguments', ...
    'B has to ba a matrix');

%% Check data in eqn structure.
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
                            assert(colE == size(B, 1), ...
                                'MESS:error_arguments', ...
                                ['number of columns of E_ differs ' ...
                                'with number of rows of B']);
                            C = eqn.EL * B;
                        case 'T' % Implement operation EL*B'.
                            assert(colE == size(B, 2), ...
                                'MESS:error_arguments', ...
                                ['number of columns of E_ differs ' ...
                                'with number of columns of B']);
                            C = eqn.EL * B';
                    end

                case 'T'
                    switch opB
                        case 'N' % Implement operation EL'*B.
                            assert(rowE == size(B, 1), ...
                                'MESS:error_arguments', ...
                                ['number of rows of E_ differs with ' ...
                                'number rows of B']);
                            C = eqn.EL' * B;
                        case 'T' % Implement operation EL'*B'.
                            assert(rowE == size(B, 2), ...
                                'MESS:error_arguments', ...
                                ['number of rows of E_ differs with ' ...
                                'number of columns of B']);
                            C = eqn.EL' * B';
                    end

            end
        case 'U'
            switch opFac
                case 'N'
                    switch opB
                        case 'N' % Implement operation EU*B.
                            assert(colE == size(B, 1), ...
                                'MESS:error_arguments', ...
                                ['number of columns of E_ differs ' ...
                                'with number of rows of B']);
                            C = eqn.EU * B;
                        case 'T' % Implement operation EU*B'.
                            assert(colE == size(B, 2), ...
                                'MESS:error_arguments', ...
                                ['number of columns of E_ differs ' ...
                                'with number of columns of B']);
                            C = eqn.EU * B';
                    end

                case 'T'
                    switch opB
                        case 'N' % Implement operation EU'*B.
                            assert(rowE == size(B, 1), ...
                                'MESS:error_arguments', ...
                                ['number of rows of E_ differs with ' ...
                                'number rows of B']);
                            C = eqn.EU' * B;
                        case 'T' % Implement operation EU'*B'.
                            assert(rowE == size(B, 2), ...
                                'MESS:error_arguments', ...
                                ['number of rows of E_ differs with ' ...
                                'number of columns of B']);
                            C = eqn.EU' * B';
                    end

            end
    end
else % Case of E_ = I_n, was set by init.
    C = B;
end
