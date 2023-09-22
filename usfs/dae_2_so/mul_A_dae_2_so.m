function C = mul_A_dae_2_so(eqn, opts, opA, B, opB)
%% function mul_A performs operation C = opA(A_)*opB(B)
%
% Input:
%   eqn     structure contains field A_
%
%   opts    struct contains parameters for the algorithm
%
%   opA     character specifies the form of opA(A_)
%           opA = 'N' performs A_*opB(B)
%           opA = 'T' performs A_'*opB(B)
%
%   B       m-x-p matrix
%
%   opB     character specifies the form of opB(B)
%           opB = 'N' performs opA(A_)*B
%           opB = 'T' performs opA(A_)*B'

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%     A = [ 0 I 0;
%           K D G';
%           G 0 0]
%
% Output:
% C = opA(A_)*opB(B)
%
%   uses size_dae_1

%% check input Parameters
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
for mat = 'EKG'
    if not(isfield(eqn, sprintf('%c_', mat))) || ...
       not(eval(sprintf('isnumeric(eqn.%c_)', mat)))
        mess_err(opts, 'error_arguments', 'field eqn.%c_ is not defined', mat);
    end
end

nv = size(eqn.M_, 1);
np = size(eqn.G_, 1);
st = 2 * nv;
n = st + np;

[rowB, colB] = size(B);

if opB == 'N'
    if n > rowB
        B = [B; zeros(np, colB)];
    elseif n < rowB
        mess_err(opts, 'error_arguments', 'B has more rows than A');
    end
else
    if n > colB
        B = [B, zeros(rowB, np)];
    elseif n < colB
        mess_err(opts, 'error_arguments', 'B has more columns than A');
    end
end

%% perform multiplication

if (opB == 'N' && (size(B, 1) == (2 * nv + np))) || (opB == 'T' && (size(B, 2) == (2 * nv + np)))
    switch opA

        case 'N'

            switch opB

                case 'N'
                    % implement operation A_*B
                    C = [B(nv + 1:2 * nv, :); ...
                         eqn.K_ * B(1:nv, :) + eqn.E_ * B(nv + 1:2 * nv, :) + eqn.G_' * B(2 * nv + 1:end, :)
                         eqn.G_ * B(1:nv, :)];
                case 'T'
                    % implement operation A_*B'
                    C = [B(:, nv + 1:2 * nv)'; ...
                         eqn.K_ * B(:, 1:nv)' + eqn.E_ * B(:, nv + 1:2 * nv)' + eqn.G_' * B(:, 2 * nv + 1:end)'
                         eqn.G_ * B(:, 1:nv)'];
            end

        case 'T'
            switch opB
                case 'N'
                    % implement operation A_'*B
                    C = [eqn.K_' * B(nv + 1:2 * nv, :) + eqn.G_' * B(2 * nv + 1:end, :); ...
                         B(1:nv, :) + eqn.E_' * B(nv + 1:2 * nv, :)
                         eqn.G_ * B(nv + 1:2 * nv, :)];
                case 'T'
                    % implement operation A_'*B'
                    C = [eqn.K_' * B(:, nv + 1:2 * nv)' + eqn.G_' * B(:, 2 * nv + 1:end)'; ...
                         B(:, 1:nv)' + eqn.E_' * B(:, nv + 1:2 * nv)'
                         eqn.G_ * B(:, nv + 1:2 * nv)'];
            end
    end

elseif (opB == 'N' && (size(B, 1) == (2 * nv))) || (opB == 'T' && (size(B, 2) == (2 * nv)))
    mess_err(opts, 'error_usage', 'mul_A_dae_2_so is only coded for shift parameter computation');
else
    mess_err(opts, 'error_arguments', 'B has wrong number of cols');
end
if opB == 'N'
    C = C(1:rowB, :);
else
    C = C(1:colB, :);
end
end
