function C = mul_E_so_1(eqn, opts, opE, B, opB)
% function C=mul_E_so_1(eqn, opts,opE,B,opB)
%
% Call help mess_usfs_so_1 to see the description of the second order
% system and its transformed first order system
%
%
% This function returns C = E*B, where matrix E given by structure
% eqn and input matrix B could be transposed.
% Matrix E is assumed to be quadratic and has a same size of 2*size(K).
%
%   Inputs:
%
%   eqn     structure containing data for matrix E (fields 'K_' and 'M_')
%   opts    structure containing parameters for the algorithm
%   opE     character specifying the shape of E
%           opE = 'N' performs E*opB(B)
%           opE = 'T' performs E'*opB(B)
%   B       m-x-p matrix
%   opB     character specifying the shape of B
%           opB = 'N' performs opE(E)*B
%           opB = 'T' performs opE(E)*B'
%
%   Output:
%
%       |-K  0|
%   C = | 0  M| *opB(B)
%
% This function does not use other so1 functions.
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
if not(isfield(eqn, 'K_')) || not(isnumeric(eqn.K_)) || ...
        not(isfield(eqn, 'M_')) || not(isnumeric(eqn.M_))
    mess_err(opts, 'error_arguments', ...
             'E consists of K and M, field eqn.K_ or eqn.M_ is not defined');
end

[rowK, colK] = size(eqn.K_);
colE = 2 * colK;

%% perform multiplication
switch opB

    % implement operation E*B
    case 'N'
        if not(colE == size(B, 1))
            mess_err(opts, 'error_arguments', ...
                     'number of columns of E differs with number of rows of B');
        end
        C = [-eqn.K_ * B(1:rowK, :)
             eqn.M_ * B(rowK + 1:end, :)];

        % implement operation E*B'
    case 'T'
        if not(colE == size(B, 2))
            mess_err(opts, 'error_arguments', ...
                     ['number of columns of E differs with number ' ...
                      'of columns of B']);
        end
        C = [-eqn.K_ * B(:, 1:colK)'; ...
             eqn.M_ * B(:, colK + 1:end)'];

end

end
