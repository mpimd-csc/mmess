function C = mul_A_dae_2(eqn, opts, opA, B, opB)
%% function mul_A performs operation C = opA(A_)*opB(B)
% Depending on the size of B either multiplication with
%     A = [A1 F;
%           G 0 ]
% or
%  A = P*A1*P'  with
%  P = I - F ( G E1\F ) \ G / E the hidden manifold projector
%  and E1 the corresponding 1,1 block in eqn.E_
%
% Input:
%   eqn     structure containing field A_ and E_
%
%   opts    struct containing parameters for the algorithm
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
%
% Output:
% C = opA(A_)*opB(B)
%

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

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
if (not(isfield(eqn, 'A_'))) || not(isnumeric(eqn.A_))
    mess_err(opts, 'error_arguments', 'field eqn.A_ is not defined');
end
if not(isfield(eqn, 'manifold_dim'))    || not(isnumeric(eqn.manifold_dim))
    mess_err(opts, 'equation_data', ...
             ['Missing or corrupted manifold_dim field detected in ' ...
              'equation structure.']);
end

n = size(eqn.A_, 1);
n_ode = eqn.manifold_dim;
one = 1:n_ode;

[rowB, colB] = size(B);

if opB == 'N'
    switch rowB
        case n
            dim = n;
        case n_ode
            dim = n_ode;
        otherwise
            mess_err(opts, 'error_arguments', ...
                     'B has wrong number of rows.');
    end
else
    switch colB
        case n
            dim = n;
        case n_ode
            dim = n_ode;
        otherwise
            mess_err(opts, 'error_arguments', ...
                     'B has wrong number of columns.');
    end
end

%% perform multiplication
if dim == n
    switch opA

        case 'N'

            switch opB

                case 'N'
                    % implement operation A_*B
                    C = eqn.A_ * B;

                case 'T'
                    % implement operation A_*B'
                    C = eqn.A_ * B';

            end

        case 'T'

            switch opB

                case 'N'
                    % implement operation A_'*B
                    C = eqn.A_' * B;

                case 'T'
                    % implement operation A_'*B'
                    C = eqn.A_' * B';

            end

    end

else

    switch opA
        case 'N'
            V = eqn.A_(one, one) * mul_Pi(eqn, opts, 'r', 'N',  B, opB);
            C = mul_Pi(eqn, opts, 'l', 'N', V, 'N');
        case 'T'
            V = eqn.A_(one, one)' * mul_Pi(eqn, opts, 'l', 'T',  B, opB);
            C = mul_Pi(eqn, opts, 'r', 'T', V, 'N');
    end
end
end
