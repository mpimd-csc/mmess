function C = mul_ApE_so_2(eqn, opts, opA, p, opE, B, opB) % function C=mul_ApE_so_2(eqn, opts,opA,p,opE,B,opB)
%
% Call help mess_usfs_so_2 to see the description of the second order
% system and its transformed first order system
%
%
% This function returns C = (A+p*E)*B, where matrices A and E given
% by structure eqn and input matrix B could be transposed.
% Matrix A is assumed to be quadratic and has a size of 2*size(K).
%
%   Inputs:
%
%   eqn     structure containing data for matrices A
%           (fields 'K_' and 'M_') and E (fields 'E_' and 'M_')
%   opts    structure containing parameters for the algorithm
%   opA     character specifying the shape of A
%           opA = 'N' performs (A + p*opE(E))*opB(B)
%           opA = 'T' performs (A' + p*opE(E))*opB(B)
%   p       scalar value
%   opE     character specifying the shape of E
%           opA = 'N' performs (opA(A) + p*E)*opB(B)
%           opA = 'T' performs (opA(A) + p*E')*opB(B)
%   B       m-x-p matrix
%   opB     character specifying the shape of B
%           opB = 'N' performs (opA(A) + p*opE(E))*B
%           opB = 'T' performs (opA(A) + p*opE(E))*B'
%
%   Output:
%
%       (|-K  0|      | E  M|)
%   C = (| 0  M| + p* | M  0|)*opB(B)
%
% This function does not use other so3 functions.
%
% ATTENTION: opA,opE are not used since matrices A and E are symmetric.

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%% check input parameters
if not(ischar(opA)) || not(ischar(opB)) || not(ischar(opE))
    mess_err(opts, 'error_arguments', 'opA, opB or opE is not a char');
end

opA = upper(opA);
opB = upper(opB);
opE = upper(opE);

if not(opA == 'N' || opA == 'T')
    mess_err(opts, 'error_arguments', 'opA is not ''N'' or ''T''');
end

if not(opB == 'N' || opB == 'T')
    mess_err(opts, 'error_arguments', 'opB is not ''N'' or ''T''');
end

if not(opE == 'N' || opE == 'T')
    mess_err(opts, 'error_arguments', 'opE is not ''N'' or ''T''');
end

if not(isnumeric(p))
    mess_err(opts, 'error_arguments', 'p is not numeric');
end

if (not(isnumeric(B))) || (not(ismatrix(B)))
    mess_err(opts, 'error_arguments', 'B has to ba a matrix');
end

%% check data in eqn structure
if not(isfield(eqn, 'haveE'))
    eqn.haveE = false;
end
if eqn.haveE
    if not(isfield(eqn, 'M_')) || not(isnumeric(eqn.M_)) || not(isfield(eqn, 'E_')) || ...
       not(isnumeric(eqn.E_)) || not(isfield(eqn, 'K_')) || not(isnumeric(eqn.K_))
        mess_err(opts, 'error_arguments', ...
                 'field eqn.M_, eqn.E_ or eqn.K_ is not defined or corrupted');
    end
else
    mess_err(opts, 'error_arguments', ['eqn.haveE has to be 1 because of ' ...
                                       'the structure of E']);
end

[rowK, colK] = size(eqn.K_);

colA = 2 * colK;

%% perform multiplication
switch opB
    % implement operation (A+p*E)*B = C
    case 'N'
        if not(colA == size(B, 1))
            mess_err(opts, 'error_arguments', ...
                     'number of columns of A differs from number of rows of B');
        end
        temp = p * eqn.M_;
        C = [(p * eqn.E_ - eqn.K_) * B(1:rowK, :) + temp * B(rowK + 1:end, :); ...
             temp * B(1:rowK, :) + eqn.M_ * B(rowK + 1:end, :)];
        % implement operation (A+p*E)*B'= C
    case 'T'
        if not(colA == size(B, 2))
            mess_err(opts, 'error_arguments', ...
                     'number of columns of A differs with number of columns of B');
        end
        temp = p * eqn.M_;
        C = [(p * eqn.E_ - eqn.K_) * B(:, 1:rowK)' + temp * B(:, rowK + 1:end)'; ...
             temp * B(:, 1:rowK)' + eqn.M_ * B(:, rowK + 1:end)'];
end
end
