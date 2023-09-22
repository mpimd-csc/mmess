function X = sol_A_dae_1_so(eqn, opts, opA, B, opB)
%% function sol_A solves opA(A) * X = opC(B) resp. performs X = opA(A) \ opB(B)
% for A as in (2) in help mess_usfs_dae1_so
%
%  X = sol_A_dae_1_so(eqn, opts, opA, B, opB)
%
% M, D, K are assumed to be quadratic.
% Input:
%
%   eqn     structure contains  data for A (E_,K_) and E (M_,K_)
%
%   opts    struct contains parameters for the algorithm
%
%   opA     character specifies the form of opA(A)
%           opA = 'N' for A
%           opA = 'T' for A'
%
%   B       n-x-p matrix
%
%   opB     character specifies the form of opB(B)
%           opB = 'N' for B
%           opB = 'T' for B'
%
% Output
%
%   X       matrix fulfills equation opA(A) * X = B
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

if not(isfield(eqn, 'K_')) || not(isnumeric(eqn.K_))
    mess_err(opts, 'equation_data', ...
             'Empty or Corrupted field K detected in equation structure.');
end

if not(isfield(eqn, 'manifold_dim')) || not(isnumeric(eqn.manifold_dim))
    mess_err(opts, 'equation_data', ...
             ['Missing or corrupted manifold_dim field detected in ' ...
              'equation structure.']);
end

n = size(eqn.K_, 1);
manifold_dim = eqn.manifold_dim;
na = n - manifold_dim;
one = 1:manifold_dim;
twob = (manifold_dim + 1):(2 * manifold_dim);

if opB == 'N'
    rows = size(B, 1);
    cols = size(B, 2);
else
    rows = size(B, 2);
    cols = size(B, 1);
end

if not(2 * manifold_dim == rows)
    mess_err(opts, 'error_arguments', ...
             ['number of rows of B differs from number of cols of A ' ...
              '(2 * manifold_dim)']);
end
%% solve

if issymmetric(eqn.K_) && issymmetric(eqn.M_)
    opA = 'N';   % let us avoid unnecessary transposition of matrices
end

switch opA

    case 'N'
        switch opB

            case 'N'
                X = eqn.K_ \ [B(one, :); zeros(na, cols)];
                X = [-X(one, :); eqn.M_(one, one) \ B(twob, :)];

            case 'T'
                X = eqn.K_ \ [B(:, one)'; zeros(na, cols)];
                X = [-X(one, :); eqn.M_(one, one) \ B(:, twob)'];
        end

    case 'T'
        switch opB

            case 'N'
                X = eqn.K_' \ [B(one, :); zeros(na, cols)];
                X = [-X(one, :); eqn.M_(one, one)' \ B(twob, :)];

            case 'T'
                X = eqn.K_' \ [B(:, one)'; zeros(na, cols)];
                X = [-X(one, :); eqn.M_(one, one)' \ B(:, twob)'];

        end

end
