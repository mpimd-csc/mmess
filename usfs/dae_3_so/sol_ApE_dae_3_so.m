function X = sol_ApE_dae_3_so(eqn, opts, opA, p, opE, B, opB)
%% function sol_ApE solves (opA(A_) + p*opE(E_))*X = opB(B)
%  resp. performs X=(opA(A_)+p*opE(E_))\opB(B)
%
%
% A_ and E_ are assumed to be quadratic.
% Input:
%
%   eqn     structure contains A_ and E_
%
%   opts    struct contains parameters for the algorithm
%
%   opA     character specifies the form of opA(A_)
%           opA = 'N' for A_
%           opA = 'T' for A_'
%
%   p       scalar Value
%
%   opE     character specifies the form of opE(E_)
%           opE = 'N' for E_
%           opE = 'T' for E_'
%
%   B       n-x-p matrix
%
%   opB     character specifies the form of opB(B)
%           opB = 'N' for B
%           opB = 'T' for B'
%
%   typeE   specifies whether E_ is Identity or not
%           typeE = 0 E_ is Identity
%           typeE = 1 E_ is not Identity
%
% Output
%
%   X       matrix fulfills equation (opA(A_)+p*opE(E_))*X = B
%

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%% check input Parameters
if not(ischar(opA)) || not(ischar(opE)) || not(ischar(opB))
    mess_err(opts, 'error_arguments', 'opA, opE or opB is not a char');
end

opA = upper(opA);
opE = upper(opE);
opB = upper(opB);

if not(opA == 'N' || opA == 'T')
    mess_err(opts, 'error_arguments', 'opA is not ''N'' or ''T''');
end

if not(opE == 'N' || opE == 'T')
    mess_err(opts, 'error_arguments', 'opE is not ''N'' or ''T''');
end

if not(opB == 'N' || opB == 'T')
    mess_err(opts, 'error_arguments', 'opB is not ''N'' or ''T''');
end

if not(isnumeric(p))
    mess_err(opts, 'error_arguments', 'p is not numeric');
end

if (not(isnumeric(B))) || (not(ismatrix(B)))
    mess_err(opts, 'error_arguments', 'B has to ba a matrix');
end

%% check data in eqn structure
for mat = 'MEKG'
    if not(isfield(eqn, sprintf('%c_', mat))) || ...
       not(eval(sprintf('isnumeric(eqn.%c_)', mat)))
        mess_err(opts, 'error_arguments', 'field eqn.%c_ is not defined', mat);
    end
end

nv = size(eqn.M_, 1);
np = size(eqn.G_, 1);

[rowB, colB] = size(B);

if opB == 'N'
    if not(rowB == 2 * nv + np)
        B = [B; zeros(2 * nv + np - rowB, colB)];
    end
else
    if not(colB == 2 * nv + np)
        B = [B, zeros(rowB, 2 * nv + np - colB)];
    end
end

switch opA

    case 'N'
        switch opE

            case 'N'

                switch opB
                    % implement solve (A_+p*E_)*X=B
                    case 'N'
                        x13 = [eqn.K_ - p * eqn.E_ - p^2 * eqn.M_, eqn.G_'; ...
                               eqn.G_, sparse(np, np)] \ ...
                              [B(nv + 1:2 * nv, :) - ...
                               eqn.E_ * B(1:nv, :) - ...
                               p * eqn.M_ * B(1:nv, :); ...
                               B(2 * nv + 1:end, :)];
                        X = [x13(1:nv, :); ...
                             B(1:nv, :) - p * x13(1:nv, :); ...
                             x13(nv + 1:end, :)];
                        % X = (A + p * E) \ B;

                        % implement solve (A_+p*E_)*X=B'
                    case 'T'
                        x13 = [eqn.K_ - p * eqn.E_ - p^2 * eqn.M_, eqn.G_'; ...
                               eqn.G_, sparse(np, np)] \ ...
                              [B(:, nv + 1:2 * nv)' - ...
                               eqn.E_ * B(:, 1:nv)' - ...
                               p * eqn.M_ * B(:, 1:nv)'; ...
                               B(:, 2 * nv + 1:end)'];
                        X = [x13(1:nv, :); B(:, 1:nv)' - ...
                             p * x13(1:nv, :); ...
                             x13(nv + 1:end, :)];
                        % X = (A + p * E) \ B';
                end

            case 'T'

                switch opB
                    % implement solve (A_+p*E_)*X=B
                    case 'N'
                        x13 = [eqn.K_ - p * eqn.E_ - p^2 * eqn.M_', eqn.G_'; ...
                               eqn.G_, sparse(np, np)] \ ...
                              [B(nv + 1:2 * nv, :) - ...
                               eqn.E_ * B(1:nv, :) - ...
                               p * eqn.M_' * B(1:nv, :); ...
                               B(2 * nv + 1:end, :)];
                        X = [x13(1:nv, :); ...
                             B(1:nv, :) - p * x13(1:nv, :); ...
                             x13(nv + 1:end, :)];
                        % X = (A + p * E') \ B;

                        % implement solve (A_+p*E_)*X=B'
                    case 'T'
                        x13 = [eqn.K_ - p * eqn.E_ - p^2 * eqn.M_', eqn.G_'; ...
                               eqn.G_, sparse(np, np)] \ ...
                              [B(:, nv + 1:2 * nv)' - ...
                               eqn.E_ * B(:, 1:nv)' - ...
                               p * eqn.M_' * B(:, 1:nv)'; ...
                               B(:, 2 * nv + 1:end)'];
                        X = [x13(1:nv, :); ...
                             B(:, 1:nv)' - p * x13(1:nv, :); ...
                             x13(nv + 1:end, :)];
                        % X = (A + p * E') \ B';
                end

        end

    case 'T'
        switch opE

            case 'N'

                switch opB
                    % implement solve (A_+p*E_)*X=B
                    case 'N'
                        x23 = [eqn.K_' - p * eqn.E_' - ...
                               p^2 * eqn.M_, eqn.G_'; ...
                               eqn.G_, sparse(np, np)] \ ...
                              [B(1:nv, :) - p * B(nv + 1:2 * nv, :); ...
                               -B(2 * nv + 1:end, :)];
                        X = [(B(1:nv, :) - ...
                              eqn.K_' * x23(1:nv, :) - ...
                              eqn.G_' * x23(nv + 1:end, :)) ./ p; ...
                             x23(1:nv, :); ...
                             x23(nv + 1:end, :)];
                        % X = (A' + p * E) \ B;

                        % implement solve (A_+p*E_)*X=B'
                    case 'T'
                        x23 = [eqn.K_' - p * eqn.E_' - ...
                               p^2 * eqn.M_, eqn.G_'; ...
                               eqn.G_, ...
                               sparse(np, np)] \ ...
                              [B(:, 1:nv)' - p * B(:, nv + 1:2 * nv)'; ...
                               -B(:, 2 * nv + 1:end)'];
                        X = [(B(:, 1:nv)' - eqn.K_' * x23(1:nv, :) - ...
                              eqn.G_' * x23(nv + 1:end, :)) ./ p; ...
                             x23(1:nv, :); ...
                             x23(nv + 1:end, :)];
                        % X = (A' + p * E) \ B';
                end

            case 'T'

                switch opB
                    % implement solve (A_+p*E_)*X=B
                    case 'N'
                        x23 = [eqn.K_' - p * eqn.E_' - p^2 * eqn.M_', ...
                               eqn.G_'; ...
                               eqn.G_, ...
                               sparse(np, np)] \ ...
                              [B(1:nv, :) - p * B(nv + 1:2 * nv, :); ...
                               -B(2 * nv + 1:end, :)];
                        X = [(B(1:nv, :) - eqn.K_' * x23(1:nv, :) - ...
                              eqn.G_' * x23(nv + 1:end, :)) ./ p; ...
                             x23(1:nv, :); ...
                             x23(nv + 1:end, :)];
                        % X = (A' + p * E') \ B;

                        % implement solve (A_+p*E_)*X=B'
                    case 'T'
                        x23 = [eqn.K_' - p * eqn.E_' - ...
                               p^2 * eqn.M_', eqn.G_'; ...
                               eqn.G_, ...
                               sparse(np, np)] \ ...
                              [B(:, 1:nv)' - p * B(:, nv + 1:2 * nv)'; ...
                               -B(:, 2 * nv + 1:end)'];
                        X = [(B(:, 1:nv)' - eqn.K_' * x23(1:nv, :) - ...
                              eqn.G_' * x23(nv + 1:end, :)) ./ p; ...
                             x23(1:nv, :); ...
                             x23(nv + 1:end, :)];
                        % X = (A' + p * E') \ B';
                end
        end

end
if opB == 'N'
    X = X(1:rowB, :);
else
    X = X(1:colB, :);
end
end
