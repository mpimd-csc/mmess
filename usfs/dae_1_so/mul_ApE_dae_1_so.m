function C = mul_ApE_dae_1_so(eqn, opts, opA, p, opE, B, opB)%#ok<INUSL>

%% function mul_A mul_ApE_dae_1_so operation C = (opA(A_)+pc*opE(E_))*opB(B)
% for A, E as in (2) in help mess_usfs_dae1_so. 
%
%  C = mul_ApE_dae_1_so(eqn, opts, opA, p, opE, B, opB)
%
% Input:
%   eqn     structure contains matrices
%
%   opts    struct contains parameters for the algorithm
%
%   opA     character specifies the form of opA(A_)
%           opA = 'N' performs (A_+pc*opE(E_))*opB(B)
%           opA = 'T' performs (A_'+pc*opE(E_))*opB(B)
%
%   p       scalar value
%
%   opE     character specifies the form of opE(E_)
%           opE = 'N' performs (opA(A_)+pc*E_)*opB(B)
%           opE = 'T' performs (opA(A_)+pc*E_')*opB(B)
%
%   B       m-x-p matrix
%
%   opB     character specifies the form of opB(B)
%           opB = 'N' performs (opA(A_)+pc*opE(E_))*B
%           opB = 'T' performs (opA(A_)+pc*opE(E_))*B'
%
% Output:
% B = (opA(A_)+pc*opE(E_))*opB(B)
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
if (not(ischar(opA)) || not(ischar(opE)) || not(ischar(opB)))
    error('MESS:error_arguments', 'opA, opE or opB is not a char');
end

opA = upper(opA); opE = upper(opE); opB = upper(opB);

if(not((opA == 'N' || opA == 'T')))
    error('MESS:error_arguments', 'opA is not ''N'' or ''T''');
end

if(not((opE == 'N' || opE == 'T')))
    error('MESS:error_arguments', 'opE is not ''N'' or ''T''');
end

if(not((opB == 'N' || opB == 'T')))
    error('MESS:error_arguments', 'opB is not ''N'' or ''T''');
end
if(not(isnumeric(p)))
    error('MESS:error_arguments','p is not numeric');
end
if (not(isnumeric(B))) || (not(ismatrix(B)))
    error('MESS:error_arguments','B has to ba a matrix');
end

%% check data in eqn structure
if (not(isfield(eqn,'K_')) || not(isnumeric(eqn.K_)))
    error('MESS:equation_data',...
        'Empty or Corrupted field K detected in equation structure.')
end
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

n = size(eqn.K_,1);
nd = eqn.nd;
one = 1 : nd;
two = (nd + 1) : n;
twob = (nd + 1) : (2 * nd);


if(opB == 'N')
    rows = size(B, 1);
else
    rows = size(B, 2);
end

if(2 * nd ~= rows)
    error('MESS:error_arguments','Rows of A differs from rows of B');
end

%% compute C = (A + p * E) * B

if issymmetric(eqn.E_) && issymmetric(eqn.M_)
    opE = 'N';   % let us avoid unnecessary transposition of matrices
end
if issymmetric(eqn.K_) && issymmetric(eqn.M_)
    opA = 'N';   % let us avoid unnecessary transposition of matrices
end

%% perfom solve operations for E ~= Identity
switch opA
    case 'N'
        switch opE
            case 'N'
                switch opB
                    case 'N'
                        C2 = eqn.M_(one, one)*...
                            (B(twob, :) + p * B(one, :));
                        C1 = -eqn.K_(one, one) * B(one,:) + ...
                            eqn.K_(one, two) * ( eqn.K_(two, two) \ ...
                            eqn.K_(two, one) * B(one,:)) + ...
                            p * (eqn.E_(one, one) * B(one, :) + ...
                            eqn.M_(one, one) * B(twob,:));
                        C = [C1; C2];
                    case 'T'
                        C2 = eqn.M_(one, one)*...
                            (B(:, twob)' + p * B(:, one)');
                        C1 = -eqn.K_(one, one) * B(:, one)' + ...
                            eqn.K_(one, two) * ( eqn.K_(two, two) \ ...
                            eqn.K_(two, one) * B(:, one)') + ...
                            p * (eqn.E_(one, one) * B(:, one)' + ...
                            eqn.M_(one, one) * B(:, twob)');
                        C = [C1; C2];
                end
            case 'T'
                switch opB
                    case 'N'
                        C2 = eqn.M_(one, one)* B(twob, :) + ...
                            p * (eqn.M_(one,one)' * B(one, :));
                        C1 = -eqn.K_(one, one) * B(one,:) + ...
                            eqn.K_(one, two) * ( eqn.K_(two, two) \ ...
                            eqn.K_(two, one) * B(one,:)) + ...
                            p * (eqn.E_(one, one)' * B(one, :) + ...
                            eqn.M_(one, one) * B(twob,:));
                        C = [C1; C2];
                    case 'T'
                        C2 = eqn.M_(one, one)* B(:, twob)' + ...
                            p * (eqn.M_(one,one)' * B(:, one)');
                        C1 = -eqn.K_(one, one) * B(:, one)' + ...
                            eqn.K_(one, two) * ( eqn.K_(two, two) \ ...
                            eqn.K_(two, one) * B(:, one)') + ...
                            p * (eqn.E_(one, one)' * B(one, :) + ...
                            eqn.M_(one, one) * B(:, twob)');
                        C = [C1; C2];
                end
        end
    case 'T'
        switch opE
            case 'N'
                switch opB
                    case 'N'
                        C2 = eqn.M_(one, one)' * B(twob, :) + ...
                            p * (eqn.M_(one, one) * B(one, :));
                        C1 = -eqn.K_(one, one)' * B(one,:) + ...
                            eqn.K_(two, one)' * ( eqn.K_(two, two)' \ ...
                            eqn.K_(one, two)' * B(one,:)) + ...
                            p * (eqn.E_(one, one) * B(one, :) + ...
                            eqn.M_(one, one) * B(twob,:));
                        C = [C1; C2];
                    case 'T'
                        C2 = eqn.M_(one, one)' * B(:, twob)' + ...
                            p * (eqn.M_(one, one) * B(:, one)');
                        C1 = -eqn.K_(one, one)' * B(:, one)' + ...
                            eqn.K_(two, one)' * (eqn.K_(two, two)' \ ...
                            eqn.K_(one, two)' * B(:, one)') + ...
                            p * (eqn.E_(one, one) * B(:, one)' + ...
                            eqn.M_(one, one) * B(:, twob)');
                        C = [C1; C2];
                end
            case 'T'
                switch opB
                    case 'N'
                        C2 = eqn.M_(one, one)'*...
                            (B(twob, :) + p * B(one, :));
                        C1 = -eqn.K_(one, one)' * B(one,:) + ...
                            eqn.K_(two, one)' * ( eqn.K_(two, two)' \ ...
                            eqn.K_(one, two)' * B(one,:)) + ...
                            p * (eqn.E_(one, one)' * B(one, :) + ...
                            eqn.M_(one, one)' * B(twob,:));
                        C = [C1; C2];
                     case 'T'
                        C2 = eqn.M_(one, one)'*...
                            (B(:, twob)' + p * B(:, one)');
                        C1 = -eqn.K_(one, one)' * B(:, one)' + ...
                            eqn.K_(two, one)' * ( eqn.K_(two, two)' \ ...
                            eqn.K_(one, two)' * B(:, one)') + ...
                            p * (eqn.E_(one, one)' * B(:, one)' + ...
                            eqn.M_(one, one)' * B(:, twob)');
                        C = [C1; C2];
                end
        end
end
