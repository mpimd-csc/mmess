function X = sol_E_so_2(eqn, opts, opE, B, opB)
% function X=sol_E_so_2(eqn, opts,opE,B,opB)
%
% Call help mess_usfs_so_2 to see the description of the second order
% system and its transformed first order system
%
%
% This function returns X= E\B, where matrix E given by structure
% eqn and input matrix B could be transposed.
%
%   Inputs:
%   eqn     structure containing data for matrix E (fields 'E_' and 'M_')
%   opts    structure containing parameters for the algorithm
%   opE     character specifying the shape of E
%           opE = 'N' solves E *X = opB(B)
%           opE = 'T' solves E'*X = opB(B)
%   B       p-x-q matrix
%   opB     character specifying the shape of B
%           opB = 'N' solves opE(E)*X = B
%           opB = 'T' solves opE(E)*X = B'
%
%   Output:
%
%                                       |E  M|
%   X       matrix fulfilling equation |M  0| *X = opB(B)
%
%   This function does not use other so3 functions.
%
% ATTENTION: opE is not used since matrix E is symmetric

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%% check input parameters
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
    mess_err(opts, 'error_arguments', 'B has to ba a matrix');
end

%% check data in eqn structure
if not(isfield(eqn, 'E_')) || not(isnumeric(eqn.E_)) || ...
        not(isfield(eqn, 'M_')) || not(isnumeric(eqn.M_))
    mess_err(opts, 'error_arguments', ...
             'E consists of D and M, field eqn.E_ or eqn.M_ is not defined');
end
if not(isfield(eqn, 'K_')) || not(isnumeric(eqn.K_))
    mess_err(opts, 'error_arguments', ...
             'Field eqn.K_ is not defined or corrupted');
end

rowK = size(eqn.K_, 1);
rowE = 2 * rowK;

%% perform solve operations
switch opB

    % implement solve E*X = B
    case 'N'

        if not(rowE == size(B, 1))
            mess_err(opts, 'error_arguments', ...
                     'number of rows of E differ with number of rows of B');
        end

        X1 = eqn.M_ \ B(rowK + 1:end, :);
        X2 = eqn.M_ \ (B(1:rowK, :) - eqn.E_ * X1);
        X  = [X1; X2];

        % implement solve E*X = B'
    case 'T'

        if not(rowE == size(B, 2))
            mess_err(opts, 'error_arguments', ...
                     'number of rows of E differs with number of columns of B');
        end

        X1 = eqn.M_ \ B(:, rowK + 1:end)';
        X2 = eqn.M_ \ (B(:, 1:rowK)' - eqn.E_ * X1);
        X  = [X1; X2];

end

end
