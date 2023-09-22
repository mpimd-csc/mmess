function C = mul_A_so_iter(eqn, opts, opA, B, opB)
% function C=mul_A_so_iter_iter(eqn,opts,opA,B,opB)
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
% This function uses another so_iter function size_so_iter_iter(eqn,
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
if not(isfield(eqn, 'M_') && isfield(eqn, 'E_') && isfield(eqn, 'K_'))
    mess_err(opts, 'error_arguments', ...
             'field eqn.M_ or eqn.E_ or eqn.K_ is not defined');
end

rowM = size_so_iter(eqn, opts);
colM = 2 * rowM;
colK = 2 * size(eqn.K_, 1);
colD = 2 * size(eqn.E_, 1);
rowB = size(B, 1);
colB = size(B, 2);
n = rowB / 2;
alpha = opts.usfs.so_iter.alpha;

%% Perform multiplication

switch opB

    % Implement multiplication A_*B
    case 'N'
        if not(colM == size(B, 1))  % M and B/2 and K should agree in size
            mess_err(opts, 'error_arguments', ...
                     ['number of columns of A differs from number ' ...
                      'of rows of B']);
        end
        if not(colK == size(B, 1))
            mess_err(opts, 'error_arguments', ...
                     ['number of columns of A differs from number ' ...
                      'of rows of B']);
        end
        if not(colD == size(B, 1))
            mess_err(opts, 'error_arguments', ...
                     ['number of columns of A differs from number ' ...
                      'of rows of B']);
        end
        B_mat = [B(1:n, :), B(n + 1:rowB, :)];
        K_mat = eqn.K_ * B_mat;
        C = [-alpha * K_mat(:, 1:colB) + K_mat(:, colB + 1:2 * colB) - alpha * eqn.E_ * B_mat(:, colB + 1:2 * colB)
             -K_mat(:, 1:colB) - eqn.E_ * B_mat(:, colB + 1:2 * colB) + alpha * eqn.M_ * B_mat(:, colB + 1:2 * colB)];

        % Implement multiplication A_*B'
    case 'T'
        if not(colM == size(B, 2))
            mess_err(opts, 'error_arguments', ...
                     ['number of columns of A differs from number ' ...
                      'of columns of B']);
        end
        if not(colK == size(B, 2))
            mess_err(opts, 'error_arguments', ...
                     ['number of columns of A differs from number ' ...
                      'of columns of B']);
        end
        if not(colD == size(B, 2))
            mess_err(opts, 'error_arguments', ...
                     ['number of columns of A differs from number ' ...
                      'of columns of B']);
        end
        B_mat = [B(:, 1:colB / 2)', B(:, 1 + colB / 2:colB)'];
        K_mat = eqn.K_ * B_mat;
        C = [-alpha * K_mat(:, 1:colB) + K_mat(:, colB + 1:2 * colB) - alpha * eqn.E_ * B_mat(:, colB + 1:2 * colB)
             -K_mat(:, 1:colB) - eqn.E_ * B_mat(:, colB + 1:2 * colB) + alpha * eqn.M_ * B_mat(:, colB + 1:2 * colB)];
end
end
