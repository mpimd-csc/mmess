function X = sol_A_dae_2(eqn, opts, opA, B, opB)
%% function sol_A solves solves opA(A_)*X = opB(B)
%
% Depending on the vertical dimension of B this solves either with
%     A = [A1 F;
%           G 0 ]
% or
%  A = P*A1*P'  with
%  P = I - F ( G E1\F ) \ G / E the hidden manifold projector
%  and E1 the corresponding 1,1 block in eqn.E_
%
% Input:
%  eqn       structure with fields A_ and E_
%  opts      struct contains parameters for the algorithm
%  opA       character specifies the form of opA(A_)
%                  opA = 'N' solves A*X=opB(B)
%                  opA = 'T' solves A^T*X=opB(B)
%
%  B         p-x-q matrix
%
%  opB       character specifies the form of opB(B)
%                  opB = 'N' solves A*X=B
%                  opB = 'T' solves A*X=B^T
%
% Output:
%  X       matrix fulfills equation opA(A)X = opB(B)
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
if not(isfield(eqn, 'A_'))
    mess_err(opts, 'error_arguments', 'field eqn.A_ is not defined');
end
if not(isfield(eqn, 'manifold_dim'))    || not(isnumeric(eqn.manifold_dim))
    mess_err(opts, 'error_arguments', ...
             ['Missing or corrupted manifold_dim field detected in ' ...
              'equation structure.']);
end

n = size(eqn.A_, 1);
n_ode = eqn.manifold_dim;

[rowB, colB] = size(B);

if opB == 'N'
    switch rowB
        case n
            dim = n;
        case n_ode
            dim = n_ode;
        otherwise
            mess_err(opts, 'error_arguments', 'B has wrong number of rows.');
    end
else
    switch colB
        case n
            dim = n;
        case n_ode
            dim = n_ode;
        otherwise
            mess_err(opts, 'error_arguments', 'B has wrong number of columns.');
    end
end

%% solve
if dim == n
    switch opA

        case 'N'
            switch opB

                % implement solve A_*X=B
                case 'N'
                    if not(n == size(B, 1))
                        mess_err(opts, 'error_arguments', ...
                                 'number of rows of A_ differs with rows of B');
                    end
                    X = eqn.A_ \ B;

                    % implement solve A_*X=B'
                case 'T'
                    if not(n == size(B, 2))
                        mess_err(opts, 'error_arguments', ...
                                 'number of rows of A_ differs with cols of B');
                    end
                    X = eqn.A_ \ B';
            end

        case 'T'
            switch opB

                % implement solve A_'*X=B
                case 'N'
                    if not(n == size(B, 1))
                        mess_err(opts, 'error_arguments', ...
                                 'number of cols of A_ differs with rows of B');
                    end
                    X = eqn.A_' \ B;

                    % implement solve A_'*X=B'
                case 'T'
                    if not(n == size(B, 2))
                        mess_err(opts, 'error_arguments', ...
                                 'number of cols of A_ differs with cols of B');
                    end
                    X = eqn.A_' \ B';
            end

    end
else
    mess_err(opts, 'error_arguments', 'A is singular in these coordinates');
end
