function X = sol_E_dae_1_so(eqn, opts, opE, B, opB)%#ok<INUSL>
%% function sol_E_dae_1_so solves opE(E)*X = opB(B), i.e., 
% it performs X=opE(E)\opB(B) with E as in (2) in help mess_usfs_dae1_so 
%
% Input:
%   eqn     structure contains data for E (here M_,K_)
%
%   opts    struct contains parameters for all algorithms used
%
%   opE     character specifies the form of opE(E)
%           opE = 'N' solves E  * X = opB(B)
%           opE = 'T' sovles E' * X = opB(B)
%
%   B       p-x-q matrix, the right hand side for the solve
%
%   opB     character specifies the form of opB(B)
%           opB = 'N' solves opE(E)*X = B
%           opB = 'T' solves opE(E)*X = B'
%
% Output
%
%   X       matrix fulfills equation opE(E)*X = opB(B)
%
%   uses no other dae_1_so function
%
% See also mess_usfs_dae_1_so

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2021 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%% check input Paramters
if (not(ischar(opE)) || not(ischar(opB)))
    error('MESS:error_arguments', 'opE or opB is not a char');
end

opE = upper(opE); opB = upper(opB);
if(not((opE == 'N' || opE == 'T')))
    error('MESS:error_arguments','opE is not ''N'' or ''T''');
end

if(not((opB == 'N' || opB == 'T')))
    error('MESS:error_arguments','opB is not ''N'' or ''T''');
end
if (not(isnumeric(B))) || (not(ismatrix(B)))
    error('MESS:error_arguments','B has to ba a matrix');
end

%% check data in eqn structure
if (not(isfield(eqn,'M_')) || not(isnumeric(eqn.M_)))
    error('MESS:equation_data',...
        'Empty or Corrupted field M detected in equation structure.')
elseif (not(isfield(eqn,'E_')) || not(isnumeric(eqn.E_)))
    error('MESS:equation_data',...
        'Empty or Corrupted field D detected in equation structure.')
end

if not(isfield(eqn, 'nd'))    || not(isnumeric(eqn.nd))
    error('MESS:nd',...
    'Missing or Corrupted nd field detected in equation structure.');
end

nd = eqn.nd;
one = 1 : nd;
twob = (nd + 1) : (2 * nd);

if(opB == 'N')
    rows = size(B, 1);
else
    rows = size(B, 2);
end

if(2 * nd ~= rows)
    error('MESS:error_arguments', ...
        'number of rows of B differs from number of cols of E ( 2 * nd)');
end

if issymmetric(eqn.E_) && issymmetric(eqn.M_)
    opE = 'N';   % let us avoid unnecessary transposition of matrices
end

%% solve
switch opE

    case 'N'
        switch opB

            case 'N'
                X1 = eqn.M_(one, one) \ B(twob, :);
                X = [X1; eqn.M_(one, one) \ ...
                    (B(one , : ) - eqn.E_(one, one) * X1)];

            case 'T'
                X1 = eqn.M_(one, one) \ B(:, twob)';
                X = [X1; eqn.M_(one, one) \ ...
                    (B(:, one)' - eqn.E_(one, one) * X1)];

        end

    case 'T'
        switch opB

            case 'N'
                X1 = eqn.M_(one, one)' \ B(twob, :);
                X = [X1; eqn.M_(one, one)' \ ...
                    (B(one , : ) - eqn.E_(one, one)' * X1)];

            case 'T'
                X1 = eqn.M_(one, one)' \ B(:, twob)';
                X = [X1; eqn.M_(one, one)' \ ...
                    (B(:, one)' - eqn.E_(one, one)' * X1)];

        end
end
