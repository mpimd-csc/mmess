function X = sol_ApE_dae_1_so(eqn, opts, opA, p, opE, C, opC)
%% function sol_ApE_so_1 solves (opA(A) + p*opE(E))*X = opC(C) respectively
%  performs X=(opA(A)+p*opE(E))\opC(C) for A, E as in (2) in
%  mess_usfs_dae1_so
%
%   X = sol_ApE_dae_1_so(eqn, opts, opA, p, opE, C, opC)
%
% Input:
%
%   eqn     structure contains  data for A (here E_, K_) and E (here M_, K_)
%
%   opts    struct contains parameters for all algorithms used
%
%   opA     character specifies the form of opA(A)
%           opA = 'N' for A
%           opA = 'T' for A'
%
%   p      scalar Value
%
%   opE     character specifies the form of opE(E)
%           opE = 'N' for E
%           opE = 'T' for E'
%
%   C       n-x-p matrix
%
%   opC     character specifies the form of opC(C)
%           opC = 'N' for C
%           opC = 'T' for C'
%
% Output
%
%   X       matrix fulfills equation (opA(A)+p*opE(E))*X = C
%
%   uses no other dae_1_so function
%
% See also mess_usfs_dae_1_so

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%% check input Parameters
if not(ischar(opA)) || not(ischar(opE)) || not(ischar(opC))
    mess_err(opts, 'error_arguments', 'opA, opE or opC is not a char');
end

opA = upper(opA);
opE = upper(opE);
opC = upper(opC);

if not(opA == 'N' || opA == 'T')
    mess_err(opts, 'error_arguments', 'opA is not ''N'' or ''T''');
end

if not(opE == 'N' || opE == 'T')
    mess_err(opts, 'error_arguments', 'opE is not ''N'' or ''T''');
end

if not(opC == 'N' || opC == 'T')
    mess_err(opts, 'error_arguments', 'opC is not ''N'' or ''T''');
end
if not(isnumeric(p))
    mess_err(opts, 'error_arguments', 'p is not numeric');
end
if (not(isnumeric(C))) || (not(ismatrix(C)))
    mess_err(opts, 'error_arguments', 'C has to ba a matrix');
end

%% check data in eqn structure
if not(isfield(eqn, 'K_')) || not(isnumeric(eqn.K_))
    mess_err(opts, 'equation_data', ...
             'Empty or Corrupted field K detected in equation structure.');
end
if not(isfield(eqn, 'M_')) || not(isnumeric(eqn.M_))
    mess_err(opts, 'equation_data', ...
             'Empty or Corrupted field M detected in equation structure.');
elseif not(isfield(eqn, 'E_')) || not(isnumeric(eqn.E_))
    mess_err(opts, 'equation_data', ...
             'Empty or Corrupted field D detected in equation structure.');
end

if not(isfield(eqn, 'manifold_dim'))    || not(isnumeric(eqn.manifold_dim))
    mess_err(opts, 'nd', ...
             ['Missing or corrupted manifold_dim field detected in ' ...
              'equation structure.']);
end
if not(isfield(eqn, 'haveE'))
    eqn.haveE = false;
end

n = size(eqn.K_, 1);
manifold_dim = eqn.manifold_dim;
one = 1:manifold_dim;
twoc = (manifold_dim + 1):(2 * manifold_dim);

if opC == 'N'
    rows = size(C, 1);
    cols = size(C, 2);
else
    rows = size(C, 2);
    cols = size(C, 1);
end

if not(2 * manifold_dim == rows)
    mess_err(opts, 'error_arguments', 'Rows of A differs from rows of C');
end

%% solve (A + p * E) * x = C
%% perform solve operations for E not the Identity
switch opA

    case 'N'
        switch opE

            case 'N'

                switch opC

                    case 'N'
                        X1 = (eqn.K_ + p * (p * eqn.M_ - eqn.E_)) \ ...
                            [p * C(twoc, :) - C(one, :); ...
                             zeros(n - manifold_dim, cols)];
                        X1 = X1(one, :);
                        X2 = eqn.M_(one, one) \ C(twoc, :) - p * X1;
                        X = [X1; X2];

                    case 'T'
                        X1 = (eqn.K_ + p * (p * eqn.M_ - eqn.E_)) \ ...
                            [p * C(:, twoc)' - C(:, one)'; ...
                             zeros(n - manifold_dim, cols)];
                        X1 = X1(one, :);
                        X2 = eqn.M_(one, one) \ C(:, twoc)' - p * X1;
                        X = [X1; X2];

                end

            case 'T'

                if not(issymmetric(eqn.M_))
                    mess_err(opts, 'notimplemented', ...
                             ['this combination of opA, opE is only ', ...
                              'available for M11 symmetric.']);
                    % this would imply eqn.M_'*eqn.M_\eqn.M_' in X1 below
                else
                    switch opC
                        case 'N'
                            X1 = (eqn.K_ + p * (p * eqn.M_ - eqn.E_')) \ ...
                                [p * C(twoc, :) - C(one, :); ...
                                 zeros(n - manifold_dim, cols)];
                            X1 = X1(one, :);
                            X2 = eqn.M_(one, one) \ C(twoc, :) - p * X1;
                            X = [X1; X2];

                        case 'T'
                            X1 = (eqn.K_ + p * (p * eqn.M_ - eqn.E_')) \ ...
                                [p * C(:, twoc)' - C(:, one)'; ...
                                 zeros(n - manifold_dim, cols)];
                            X1 = X1(one, :);
                            X2 = eqn.M_(one, one) \ C(:, twoc)' - p * X1;
                            X = [X1; X2];

                    end

                end

        end

    case 'T'
        switch opE

            case 'N'

                if not(issymmetric(eqn.M_))
                    mess_err(opts, 'notimplemented', ...
                             ['this combination of opA, opE is only ', ...
                              'available for M11 symmetric.']);
                    % this would imply eqn.M_*eqn.M_'\eqn.M_ in X1 below

                else

                    switch opC

                        case 'N'
                            X1 = (eqn.K_' + p * (p * eqn.M_ - eqn.E_)) \ ...
                                [p * C(twoc, :) - C(one, :); ...
                                 zeros(n - manifold_dim, cols)];
                            X1 = X1(one, :);
                            X2 = eqn.M_(one, one) \ C(twoc, :) - p * X1;
                            X = [X1; X2];

                        case 'T'
                            X1 = (eqn.K_' + p * (p * eqn.M_ - eqn.E_)) \ ...
                                [p * C(:, twoc)' - C(:, one)'; ...
                                 zeros(n - manifold_dim, cols)];
                            X1 = X1(one, :);
                            X2 = eqn.M_(one, one) \ C(:, twoc)' - p * X1;
                            X = [X1; X2];

                    end

                end

            case 'T'

                switch opC

                    case 'N'
                        X1 = (eqn.K_' + p * (p * eqn.M_' - eqn.E_')) \ ...
                            [p * C(twoc, :) - C(one, :); ...
                             zeros(n - manifold_dim, cols)];
                        X1 = X1(one, :);
                        X2 = eqn.M_(one, one)' \ C(twoc, :) - p * X1;
                        X = [X1; X2];

                    case 'T'
                        X1 = (eqn.K_' + p * (p * eqn.M_' - eqn.E_')) \ ...
                            [p * C(:, twoc)' - C(:, one)'; ...
                             zeros(n - manifold_dim, cols)];
                        X1 = X1(one, :);
                        X2 = eqn.M_(one, one)' \ C(:, twoc)' - p * X1;
                        X = [X1; X2];

                end

        end

end
