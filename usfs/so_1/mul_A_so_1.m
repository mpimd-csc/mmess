function C = mul_A_so_1(eqn, opts, opA, B, opB)
% function C = mul_A_so_1(eqn, opts, opA, B, opB)
%
% Call help mess_usfs_so_1 to see the description of the second order
% system and its transformed first order system
%
%
% This function returns C = A*B, where matrix A given by structure
% eqn and input matrix B could be transposed.
% Matrix A is assumed to be quadratic and has a size of 2 * size(K).
%
%   Inputs:
%
%   eqn     structure containing  data for matrix A (fields 'K_' and 'E_')
%   opts    structure containing parameters for the algorithm
%   opA     character specifying the shape of A
%           opA = 'N' performs A*opB(B)
%           opA = 'T' performs A'*opB(B)
%   B       m-x-p matrix
%   opB     character specifying the shape of B
%           opB = 'N' performs opA(A)*B
%           opB = 'T' performs opA(A)*B'
%
%   Output:
%
%       | 0 -K|
%   C = |-K -E| *opB(B)
%
% This function does not use other so1 functions.
%
% ATTENTION: opA is not used since matrix A is symmetric

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%% check input parameters
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
    mess_err(opts, 'error_arguments', 'B has to ba a matrix');
end

%% check data in eqn structure
if  not(isfield(eqn, 'K_')) || not(isnumeric(eqn.K_)) || ...
    not(isfield(eqn, 'E_')) || not(isnumeric(eqn.E_))
    mess_err(opts, 'error_arguments', ...
             ['A consists of K and E, field eqn.K_ or eqn.E_ is not ' ...
              'defined or corrupted']);
end

[rowK, colK] = size(eqn.K_);
colA = 2 * colK;

%% perform multiplication
switch opB

    % implement operation A*B
    case 'N'
        if not(colA == size(B, 1))
            mess_err(opts, 'error_arguments', ...
                     'number of columns of A differs with number of rows of B');
        end
        C = [-eqn.K_ * B(rowK + 1:end, :)
             -eqn.K_ * B(1:rowK, :) - eqn.E_ * B(rowK + 1:end, :)];

        % implement operation A*B'
    case 'T'
        if not(colA == size(B, 2))
            mess_err(opts, 'error_arguments', ...
                     ['number of columns of A differs with number ' ...
                      'of columns of B']);
        end
        C = [-eqn.K_ * B(:, colK + 1:end)'; ...
             -eqn.K_ * B(:, 1:colK)' - eqn.E_ * B(:, colK + 1:end)'];

end

end
