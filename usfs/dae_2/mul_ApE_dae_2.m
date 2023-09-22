function C = mul_ApE_dae_2(eqn, opts, opA, p, opE, B, opB)

%% function mul_ApE_default performs operation C = (opA(A_)+p*opE(E_))*opB(B)
%
% Input:
%   eqn     structure contains A_
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
% Output:
% C = (opA(A_)+p*opE(E_))*opB(B)
%
%   uses size_dae_2

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%% check input Parameters
if not(ischar(opA)) || not(ischar(opB)) || not(ischar(opE))
    mess_err(opts, 'error_arguments', 'opA, opE or opB is not a char');
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
if (not(isfield(eqn, 'A_'))) || not(isnumeric(eqn.A_))
    mess_err(opts, 'error_arguments', 'field eqn.A_ is not defined');
end

if (not(isfield(eqn, 'E_'))) || not(isnumeric(eqn.E_))
    mess_err(opts, 'error_arguments', 'field eqn.E_ is not defined');
end
n = size(eqn.A_, 1);
n_ode = eqn.manifold_dim;

[rowB, colB] = size(B);

if opB == 'N'
    if n > rowB
        B = [B; zeros(n - n_ode, colB)];
    elseif n < rowB
        mess_err(opts, 'error_arguments', 'B has more rows than A');
    end
else
    if n > colB
        B = [B, zeros(rowB, n - n_ode)];
    elseif n < colB
        mess_err(opts, 'error_arguments', 'B has more columns than A');
    end
end

%% perform multiplication
switch opA

    case 'N'

        switch opB

            case 'N'
                % implement operation (A_+p*E_)*B
                C = (eqn.A_ + p * eqn.E_) * B;

            case 'T'
                % implement operation (A_+p*E_)*B'
                C = (eqn.A_ + p * eqn.E_) * B';
        end

    case 'T'

        switch opB

            case 'N'
                % implement operation (A_+p*E_)'*B
                C = (eqn.A_ + p * eqn.E_)' * B;

            case 'T'
                % implement operatio (A_+p*E_)'*B'
                C = (eqn.A_ + p * eqn.E_)' * B';
        end

end
end
