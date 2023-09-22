function X = sol_ApE_dae_2(eqn, opts, opA, p, opE, B, opB)

%% function sol_ApE solves (opA(A_) + p*opE(E_))*X = opB(B) resp. performs X=(opA(A_)+p*opE(E_))\opB(B)
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
if not(isfield(eqn, 'A_')) || not(isnumeric(eqn.A_))
    mess_err(opts, 'equation_data', ...
             'Empty or Corrupted field A detected in equation structure.');
end
if not(isfield(eqn, 'E_')) || not(isnumeric(eqn.E_))
    mess_err(opts, 'equation_data', ...
             'Empty or Corrupted field E detected in equation structure.');
end
if not(isfield(eqn, 'manifold_dim'))    || not(isnumeric(eqn.manifold_dim))
    mess_err(opts, 'error_arguments', ...
             ['Missing or corrupted manifold_dim field detected in ' ...
              'equation structure.']);
end
if not(isfield(eqn, 'haveE'))
    eqn.haveE = false;
end
n = size(eqn.A_, 1);
n_ode = eqn.manifold_dim;

[rowB, colB] = size(B);

if opB == 'N'
    if not(rowB == n_ode)
        mess_err(opts, 'error_arguments', 'B has not same number of rows as A');
    end
    B = [B; zeros(n - n_ode, colB)];
else
    if not(colB == n_ode)
        mess_err(opts, 'error_arguments', 'B has not same number of rows as A');
    end
    B = [B, zeros(rowB, n - n_ode)];
end

%% perform solve operations for not(E_ == Identity)
if eqn.haveE
    switch opA

        case 'N'
            switch opE

                case 'N'

                    switch opB

                        % implement solve (A_+p*E_)*X=B
                        case 'N'
                            X = (eqn.A_ + p * eqn.E_) \ B;

                            % implement solve (A_+p*E_)*X=B'
                        case 'T'
                            X = (eqn.A_ + p * eqn.E_) \ B';

                    end

                case 'T'

                    switch opB

                        % implement solve (A_+p*E_')*X=B
                        case 'N'
                            X = (eqn.A_ + p * eqn.E_') \ B;

                            % implement solve (A_+p*E_')*X=B'
                        case 'T'
                            X = (eqn.A_ + p * eqn.E_') \ B';

                    end

            end

        case 'T'
            switch opE

                case 'N'

                    switch opB

                        % implement solve (A_'+p*E_)*X=B
                        case 'N'
                            X = (eqn.A_' + p * eqn.E_) \ B;

                            % implement solve (A_'+p*E_)*X=B'
                        case 'T'
                            X = (eqn.A_' + p * eqn.E_) \ B';

                    end

                case 'T'

                    switch opB

                        % implement solve (A_'+p*E_')*X=B
                        case 'N'
                            X = (eqn.A_' + p * eqn.E_') \ B;

                            % implement solve (A_'+p*E_')*X=B'
                        case 'T'
                            X = (eqn.A_' + p * eqn.E_') \ B';

                    end
            end

    end
elseif not(eqn.haveE)
    %% perform solve operations for E_ = Identity
    switch opA

        case 'N'

            switch opB

                % implement solve (A_+p*E_)*X=B
                case 'N'
                    X = (eqn.A_ + p * eqn.E_) \ B;

                    % implement solve (A_+p*E_)*X=B'
                case 'T'
                    X = (eqn.A_ + p * eqn.E_) \ B';

            end

        case 'T'

            switch opB

                % implement solve (A_'+p*E_)*X=B
                case 'N'
                    X = (eqn.A_' + p * eqn.E_) \ B;

                    % implement solve (A_'+p*E_)*X=B'
                case 'T'
                    X = (eqn.A_' + p * eqn.E_) \ B';

            end

    end
end
X = X(1:n_ode, :);
end
