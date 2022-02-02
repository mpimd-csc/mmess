function C=mul_E_dae_1_so(eqn, opts, opE, B, opB)%#ok<INUSL>

%% function mul_A_so_1 performs operation C = opE(E)*opB(B)
% for E as in (2) in help mess_usfs_dae1_so 
%
%  C = mul_E_dae_1_so(eqn, opts, opE, B, opB)
%
% Input:
%   eqn     structure contains field E
%
%   opts    struct contains parameters for the algorithm
%
%   opE     character specifies the form of opE(E)
%           opE = 'N' performs E*opB(B)
%           opE = 'T' performs E'*opB(B)
%
%   B       m-x-p matrix
%
%   opB     character specifies the form of opB(B)
%           opB = 'N' performs opE(E)*B
%           opB = 'T' performs opE(E)*B'
%
% Output:
% C = opE(E)*opB(B)
%
%   uses no other dae_1_so function
%
% See also mess_usfs_dae_1_so

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%


%% check input Parameters
if (not(ischar(opE)) || not(ischar(opB)))
    error('MESS:error_arguments', 'opE or opB is not a char');
end

opE = upper(opE); opB = upper(opB);
if(not((opE=='N' || opE=='T')))
    error('MESS:error_arguments','opE is not ''N'' or ''T''');
end

if(not((opB=='N' || opB=='T')))
    error('MESS:error_arguments','opB is not ''N'' or ''T''');
end
if (not(isnumeric(B))) || (not(ismatrix(B)))
    error('MESS:error_arguments','B has to ba a matrix');
end

%% check data in eqn structure
if (not(isfield(eqn,'M_')) || not(isnumeric(eqn.M_)))
    error('MESS:equation_data', ...
        'Empty or Corrupted field M detected in equation structure.')
elseif (not(isfield(eqn,'E_')) || not(isnumeric(eqn.E_)))
    error('MESS:equation_data', ...
        'Empty or Corrupted field D detected in equation structure.')
end

if not(isfield(eqn, 'nd'))    || not(isnumeric(eqn.nd))
    error('MESS:nd', ...
        'Missing or Corrupted nd field detected in equation structure.');
end

nd = eqn.nd;
one = 1 : nd;
twob = (nd + 1) : (2 * nd);

if(opB == 'N')
    [rows, ~] = size(B);
else
    [~, rows] = size(B);
end

if(2 * nd ~= rows)
    error('MESS:error_arguments', ...
        'number of rows of B differs from number of cols of E ( 2 * nd)');
end

if issymmetric(eqn.E_) && issymmetric(eqn.M_)
    opE = 'N';   % let us avoid unnecessary transposition of matrices
end

%% perform multiplication
switch opE

    case 'N'
        switch opB
            case 'N'
                C = [eqn.E_(one,one) * B(one,:) ...
                    + eqn.M_(one,one) * B(twob, :);
                    eqn.M_(one, one) * B(one, :)];

            case 'T'
                C = [eqn.E_(one,one) * B(:, one)' ...
                    + eqn.M_(one,one) * B(:, twob)';
                    eqn.M_(one, one) * B(:, one)'];

        end

    case 'T'
        switch opB

            case 'N'
                C = [eqn.E_(one,one)' * B(one,:) ...
                    + eqn.M_(one,one)' * B(twob, :);
                    eqn.M_(one, one)' * B(one, :)];


            case 'T'
                C = [eqn.E_(one,one)' * B(:, one)' ...
                    + eqn.M_(one,one)' * B(:, twob)';
                    eqn.M_(one, one)' * B(:, one)'];

        end

end
