function X = sol_A_so_2(eqn, opts, opA, B, opB)
% function X=sol_A_so_2(eqn, opts,opA,B,opB)
%
% Call help mess_usfs_so_2 to see the description of the second order
% system and its transformed first order system
%
%
% This function returns X = A\B, where matrix A given by structure eqn and input matrix B could be transposed.
% Matrix A is assumed to be quadratic and has a size of 2*size(K).
%
%   Inputs:
%
%   eqn     structure containing data for matrix A (fields 'K_' and 'M_')
%   opts    structure containing parameters for the algorithm
%   opA     character specifying the shape of A
%           opA = 'N' solves A *X = opB(B)
%           opA = 'T' solves A'*X = opB(B)
%   B       p-x-q matrix
%   opB     character specifying the shape of B
%           opB = 'N' solves opA(A)*X = B
%           opB = 'T' solves opA(A)*X = B'
%
%   Output:
%                                       |-K  0|
%   X       matrix fulfilling equation | 0  M| *X= opB(B)
%
%   This function does not use other so3 functions.
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
if not(isfield(eqn, 'K_')) || not(isnumeric(eqn.K_)) || ...
        not(isfield(eqn, 'M_'))  || not(isnumeric(eqn.M_))
    mess_err(opts, 'error_arguments', ...
             'A consists of K and M, field eqn.K_ or eqn.M_ is not defined');
end

rowK = size(eqn.K_, 1);
rowA = 2 * rowK;

%% perform solve operations
switch opB

    % implement solve A*X = B
    case 'N'
        if not(rowA == size(B, 1))
            mess_err(opts, 'error_arguments', ...
                     'number of rows of A differs with number of rows of B');
        end

        X1 =  eqn.K_ \ B(1:rowK, :);
        X2 =  eqn.M_ \ B(rowK + 1:end, :);
        X  =  [-X1; X2];

        % implement solve A*X = B'
    case 'T'
        if not(rowA == size(B, 2))
            mess_err(opts, 'error_arguments', ...
                     'number of rows of A differs with number of columns of B');
        end

        X1 =  eqn.K_ \ B(:, 1:rowK)';
        X2 =  eqn.M_ \ B(:, rowK + 1:end)';
        X  =  [-X1; X2];
end

end
