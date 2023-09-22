function C = mul_A_dae_1_so(eqn, opts, opA, B, opB)
%% function mul_A performs operation C = opA(A_) * opB(B)
% for A as in (2) in help mess_usfs_dae1_so
%
%  C = mul_A_dae_1_so(eqn, opts, opA, B, opB)
%
% Input:
%   eqn     structure contains field K_, E_, M_
%
%   opts    struct contains parameters for the algorithm
% s
%   opA     character specifies the form of opA(A_)
%           opA = 'N' performs A_*opB(B)
%           opA = 'T' performs A_'*opB(B)
%
%   B       2*eqn.np times p matrix
%
%   opB     character specifies the form of opB(B)
%           opB = 'N' performs opA(A_)*B
%           opB = 'T' performs opA(A_)*B'
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
    mess_err(opts, 'error_arguments', 'B has to be a matrix');
end

%% check data in eqn structure
if not(isfield(eqn, 'K_')) || not(isnumeric(eqn.K_))
    mess_err(opts, 'equation_data', ...
             'Empty or corrupted field K detected in equation structure.');
end

if not(isfield(eqn, 'manifold_dim')) || not(isnumeric(eqn.manifold_dim))
    mess_err(opts, 'equation_data', ...
             ['Missing or corrupted manifold_dim field detected in ' ...
              'equation structure.']);
end

n = size(eqn.K_, 1);
manifold_dim = eqn.manifold_dim;
one = 1:manifold_dim;
two = (manifold_dim + 1):n;
twob = (manifold_dim + 1):(2 * manifold_dim);

if opB == 'N'
    nrows = size(B, 1);
else
    nrows = size(B, 2);
end

if not(2 * manifold_dim == nrows)
    mess_err(opts, 'error_arguments', ...
             ['number of rows of B differs from number of cols of A' ...
              '(2 * manifold_dim)']);
end

if issymmetric(eqn.K_) && issymmetric(eqn.M_)
    opA = 'N';   % let us avoid unnecessary transposition of matrices
end

%% perform multiplication
switch opA

    case 'N'
        switch opB

            case 'N'
                C = [-eqn.K_(one, one) * B(one, :) + ...
                     eqn.K_(one, two) * (eqn.K_(two, two) \ ...
                                         (eqn.K_(two, one) * B(one, :)))
                     eqn.M_(one, one) * B(twob, :)];

            case 'T'
                C = [-eqn.K_(one, one) * B(:, one)' + ...
                     eqn.K_(one, two) * (eqn.K_(two, two) \ ...
                                         (eqn.K_(two, one) * B(:, one)'))
                     eqn.M_(one, one) * B(:, twob)'];
        end

    case 'T'
        switch opB

            case 'N'
                C = [-eqn.K_(one, one)' * B(one, :) + ...
                     eqn.K_(two, one)' * (eqn.K_(two, two)' \ ...
                                          (eqn.K_(one, two)' * B(one, :)))
                     eqn.M_(one, one)' * B(twob, :)];

            case 'T'
                C = [-eqn.K_(one, one)' * B(:, one)' + ...
                     eqn.K_(two, one)' * (eqn.K_(two, two)' \ ...
                                          (eqn.K_(one, two)' * B(:, one)'))
                     eqn.M_(one, one)' * B(:, twob)'];
        end

end
