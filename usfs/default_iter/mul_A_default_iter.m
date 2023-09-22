function C = mul_A_default_iter(eqn, opts, opA, B, opB)
% function C=mul_A_default_iter(eqn,opts,opA,B,opB)
%
% This function returns C = A_*B, where matrix A_ given by
% structure eqn and input matrix B could be transposed.
% Matrix A_ is assumed to be quadratic.
%
%   Inputs:
%
%   eqn     structure containing field 'A_'
%   opts    structure containing parameters for the algorithm
%   opA     character specifying the shape of A_
%           opA = 'N' performs A_*opB(B)
%           opA = 'T' performs A_'*opB(B)
%   B       m-x-p matrix
%   opB     character specifying the shape of B
%           opB = 'N' performs opA(A_)*B
%           opB = 'T' performs opA(A_)*B'
%
%   Output:
%
%   C = opA(A_)*opB(B)
%
% This function uses another default function size_default_iter(eqn,
% opts) to obtain the number of rows of matrix A_ in structure eqn,
% that should be equal to the number of rows of matrix E_.
%

% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%
%% Check input parameters
if not(ischar(opA)) || not(ischar(opB))
    mess_err(opts, 'error_arguments', 'opA or opB is not a char');
end

opA = upper(opA);
opB = upper(opB);
if not(opA == 'N' || opA == 'T')
    mess_err(opts, 'error_arguments', 'opA is not ''N'' or ''T''');
end

if not(opB == 'N' || opB == 'T')
    mess_err(opts, 'error_arguments', 'opB is not ''N'' or ''T''');
end

if (not(isnumeric(B))) || (not(ismatrix(B)))
    mess_err(opts, 'error_arguments', 'B has to be a matrix');
end

%% Check data in eqn structure
if not(isfield(eqn, 'A_'))
    mess_err(opts, 'error_arguments', 'field eqn.A_ is not defined');
end

rowA = size_default_iter(eqn, opts);
colA = rowA;

%% Perform multiplication
switch opA

    case 'N'
        switch opB

            % Implement multiplication A_*B
            case 'N'
                if not(colA == size(B, 1))
                    mess_err(opts, 'error_arguments', ...
                             ['number of columns of A_ ' ...
                              'differs from number of rows of B']);
                end
                C = eqn.A_ * B;

                % Implement multiplication A_*B'
            case 'T'
                if not(colA == size(B, 2))
                    mess_err(opts, 'error_arguments', ...
                             ['number of columns of A_ ' ...
                              'differs from number of columns of B']);
                end
                C = eqn.A_ * B';
        end

    case 'T'
        switch opB

            % Implement multiplication A_'*B
            case 'N'
                if not(rowA == size(B, 1))
                    mess_err(opts, 'error_arguments', ...
                             ['number of rows of A_ ' ...
                              'differs from number rows of B']);
                end
                C = eqn.A_' * B;

                % Implement multiplication A_'*B'
            case 'T'
                if not(rowA == size(B, 2))
                    mess_err(opts, 'error_arguments', ...
                             ['number of rows of A_ differs from ' ...
                              'number of columns of B']);
                end
                C = eqn.A_' * B';
        end

end
end
