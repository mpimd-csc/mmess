function C = mul_E_so_iter(eqn, opts, opE, B, opB)

% function C=mul_E_so_iter(eqn, opts,opE,B,opB)
%
% This function returns C = E_*B, where matrix E_ given by structure
% eqn and input matrix B could be transposed. Matrix E_ is assumed
% to be quadratic and has the same size as A_ in structure eqn.
%
%   Inputs:
%
%   eqn     structure containing field 'E_'
%   opts    structure containing parameters for the algorithm
%   opE     character specifying the shape of E_
%           opE = 'N' performs E_*opB(B)
%           opE = 'T' performs E_'*opB(B)
%   B       m-x-p matrix
%   opB     character specifying the shape of B
%           opB = 'N' performs opE(E_)*B
%           opB = 'T' performs opE(E_)*B'
%
%   Output:
%
%   C = opE(E_)*opB(B)
%
% This function uses another so_iter function
% size_so_iter_iter(eqn,opts) to obtain the number of rows of matrix A_
% in structure eqn, that should be equal to the number of rows of
% the matrix E_.

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%% Check input parameters
if not(ischar(opE)) || not(ischar(opB))
    mess_err(opts, 'error_arguments', 'opE or opB is not a char');
end

opE = upper(opE);
opB = upper(opB);
if not(opE == 'N' || opE == 'T')
    mess_err(opts, 'error_arguments', 'opE is not ''N'' or ''T''');
end

if not(opB == 'N' || opB == 'T')
    mess_err(opts, 'error_arguments', 'opB is not ''N'' or ''T''');
end

if (not(isnumeric(B))) || (not(ismatrix(B)))
    mess_err(opts, 'error_arguments', 'B has to be a matrix');
end

%% Check data in eqn structure
if not(isfield(eqn, 'M_') && isfield(eqn, 'K_'))
    mess_err(opts, 'error_arguments', 'field eqn.M_ or eqn.K_is not defined');
end

if not(mod(size(B, 1), 2) == 0)
    mess_err(opts, 'error_arguments', ...
             'matrix B must have an even number of rows');
end

rowM = size_so_iter(eqn, opts);
colM = 2 * rowM;
colK = 2 * size(eqn.K_, 1);
rowB = size(B, 1);
colB = size(B, 2);
n = rowB / 2;
alpha = opts.usfs.so_iter.alpha;

%% Perform multiplication
% Since the implicit E is symmetric there is no need to do this switch

switch opB

    % Implement multiplication E_*B
    case 'N'
        if not(colM == size(B, 1))  % M and B/2 and K should agree in size
            mess_err(opts, 'error_arguments', ...
                     ['number of columns of M_ differs from number ' ...
                      'of rows of B']);
        end
        if not(colK == size(B, 1))
            mess_err(opts, 'error_arguments', ...
                     ['number of columns of K_ differs from number ' ...
                      'of rows of B']);
        end
        % B_mat = reshape(B,[],2*size(B,2)); % reshape places the columns
        % first the top and right next the button.
        B_mat = [B(1:n, :), B(n + 1:rowB, :)];
        M_mat = eqn.M_ * B_mat;
        C = [eqn.K_ * B_mat(:, 1:colB) + alpha * M_mat(:, colB + 1:2 * colB)
             alpha * M_mat(:, 1:colB) + M_mat(:, colB + 1:2 * colB)];
        % Implement multiplication E_*B'
    case 'T'
        if not(colM == size(B, 2))
            mess_err(opts, 'error_arguments', ...
                     ['number of columns of M_ differs from number ' ...
                      'of columns of B']);
        end
        if not(colK == size(B, 2))
            mess_err(opts, 'error_arguments', ...
                     ['number of columns of K_ differs from number ' ...
                      'of columns of B']);
        end
        B_mat = [B(:, 1:colB / 2)', B(:, 1 + colB / 2:colB)'];
        M_mat = eqn.M_ * B_mat;
        C = [eqn.K_ * B_mat(:, 1:colB) + alpha * M_mat(:, colB + 1:2 * colB)
             alpha * M_mat(:, 1:colB) + M_mat(:, colB + 1:2 * colB)];
end

end
