function X = sol_A_default_iter(eqn, opts, opA, B, opB)
% function X=sol_A_default_iter(eqn, opts,opA,B,opB)
%
% This function returns X = A_\B, where matrix A_ given by
% structure eqn and input matrix B could be transposed. Matrix A_
% is assumed to be quadratic.
%
%    Inputs:
%
%    eqn       structure containing field 'A_'
%    opts      structure containing parameters for the algorithm
%    opA       character specifying the shape of A_
%                  opA = 'N' solves  A_*X = opB(B)
%                  opA = 'T' solves  A_'*X = opB(B)
%    B         p-x-q matrix
%    opB       character specifying the shape of B
%                  opB = 'N' solves  opA(A_)*X = B
%                  opB = 'T' solves  opA(A_)*X = B'
%
%    Output:
%
%    X         matrix fulfilling equation  opA(A_)*X = opB(B)
%
% This function uses another default function size_default_iter(eqn,
% opts) to obtain the number of rows of matrix A_ in structure eqn,
% that should be equal to the number of rows of matrix E_.

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%
%% Check input parameters
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

% Initial guess for vector X

if not(isfield(opts.usfs.default_iter, 'X0_A'))
    opts.usfs.default_iter.X0_A = [];
end

%% Check data in eqn structure
switch lower(opts.usfs.default_iter.method_A)
    case {'minres', 'pcg', 'symmlq', 'pcr'}
        if not(opts.usfs.default_iter.A_is_sym)
            mess_err(opts, 'error_arguments', ...
                     'field eqn.A_ is not symmetric.');
        end
end

n = size_default_iter(eqn, opts);

%% Preallocate solution
if opB == 'N'
    X = zeros(size(B));
    flags = zeros(1, size(B, 2));
else
    X = zeros(size(B'));
    flags = zeros(1, size(B, 1));
end

%% Create anonymous functions
% To call multiplication with A respecting opA
switch lower(opts.usfs.default_iter.method_A)
    case {'bicg', 'lsqr', 'qmr'}

        mul_A = @(X, flag) flagged_mul_A(X, flag, eqn, opts, opA);
    case 'pcg'

        mul_A = @(X) -mul_A_default_iter(eqn, opts, opA, X, 'N');
    otherwise

        mul_A = @(X) mul_A_default_iter(eqn, opts, opA, X, 'N');
end

% For calling the actual iterative solver

solver = eval(sprintf('@%s', lower(opts.usfs.default_iter.method_A)));

%% Perform solve operations
switch opB

    case 'N'
        if not(n == size(B, 1))
            mess_err(opts, 'error_arguments', ...
                     ['Number of rows of A_ differs from number ' ...
                      'of rows of B']);
        end

        switch lower(opts.usfs.default_iter.method_A)

            case 'pcg'
                for i = 1:size(B, 2)
                    [x, flags(i)] = ...
                        solver(mul_A, B(:, i), ...
                               opts.usfs.default_iter.res_tol, ...
                               opts.usfs.default_iter.max_iter, ...
                               opts.usfs.default_iter.PA_L, ...
                               opts.usfs.default_iter.PA_R, ...
                               opts.usfs.default_iter.X0_A);
                    X(:, i) = -x;
                end
            case 'gmres'
                for i = 1:size(B, 2)
                    [X(:, i), flags(i)] = ...
                        solver(mul_A, B(:, i), ...
                               opts.usfs.default_iter.restIter, ...
                               opts.usfs.default_iter.res_tol, ...
                               opts.usfs.default_iter.max_iter, ...
                               opts.usfs.default_iter.PA_L, ...
                               opts.usfs.default_iter.PA_R, ...
                               opts.usfs.default_iter.X0_A);
                end
            case 'pcr'
                for i = 1:size(B, 2)
                    [X(:, i), flags(i)] = ...
                        solver(mul_A, B(:, i), ...
                               opts.usfs.default_iter.res_tol, ...
                               opts.usfs.default_iter.max_iter, ...
                               @(X) mfun(X, opts), ...
                               opts.usfs.default_iter.X0_A);
                end
            otherwise
                for i = 1:size(B, 2)
                    [X(:, i), flags(i)] = ...
                        solver(mul_A, B(:, i), ...
                               opts.usfs.default_iter.res_tol, ...
                               opts.usfs.default_iter.max_iter, ...
                               opts.usfs.default_iter.PA_L, ...
                               opts.usfs.default_iter.PA_R, ...
                               opts.usfs.default_iter.X0_A);
                end

        end

    case 'T'
        if not(n == size(B, 2))
            mess_err(opts, 'error_arguments', ...
                     ['Number of rows of A_ differs from number ' ...
                      'of columns of B']);
        end

        switch lower(opts.usfs.default_iter.method_A)

            case 'pcg'
                for i = 1:size(B, 1)
                    [x, flags(i)] = ...
                        solver(mul_A, B(i, :)', ...
                               opts.usfs.default_iter.res_tol, ...
                               opts.usfs.default_iter.max_iter, ...
                               opts.usfs.default_iter.PA_L, ...
                               opts.usfs.default_iter.PA_R, ...
                               opts.usfs.default_iter.X0_A);
                    X(:, i) = -x;
                end
            case 'gmres'
                for i = 1:size(B, 1)
                    [X(:, i), flags(i)] = ...
                        solver(mul_A, B(i, :)', ...
                               opts.usfs.default_iter.restIter, ...
                               opts.usfs.default_iter.res_tol, ...
                               opts.usfs.default_iter.max_iter, ...
                               opts.usfs.default_iter.PA_L, ...
                               opts.usfs.default_iter.PA_R, ...
                               opts.usfs.default_iter.X0_A);
                end
            case 'pcr'
                for i = 1:size(B, 1)
                    [X(:, i), flags(i)] = ...
                        solver(mul_A, B(i, :)', ...
                               opts.usfs.default_iter.res_tol, ...
                               opts.usfs.default_iter.max_iter, ...
                               @(X) mfun(X, opts), ...
                               opts.usfs.default_iter.X0_A);
                end
            otherwise
                for i = 1:size(B, 1)
                    [X(:, i), flags(i)] = ...
                        solver(mul_A, B(i, :)', ...
                               opts.usfs.default_iter.res_tol, ...
                               opts.usfs.default_iter.max_iter, ...
                               opts.usfs.default_iter.PA_L, ...
                               opts.usfs.default_iter.PA_R, ...
                               opts.usfs.default_iter.X0_A);
                end

        end

end
if any(flags)
    mess_warn(opts, 'usfs_iter', ...
              [lower(opts.usfs.default_iter.method_A) ...
               ' did not converge as desired']);
    mess_fprintf(opts, ...
                 ['These are the right hand side indices and ', ...
                  'corresponding non-zero termination flags ', ...
                  'we encountered:\n']);
    idx = find(not(flags == 0));
    for iidx = 1:length(idx)
        mess_fprintf(opts, '%d %d\n', idx(iidx), flags(idx(iidx)));
    end
end
end

function Y = flagged_mul_A(X, flag, eqn, opts, opA)
% function Y = flagged_mul_A(X, flag, eqn, opts, opA)
% This is a function handle that accepts the vector input X and the matrix
% A_ given by the eqn structure, and returns the matrix vector product
% A_*X. The input 'flag' defines whether A_ should be transposed or not.

switch lower(flag)
    case 'notransp'
        my_opA = opA;
    case 'transp'
        if strcmp(opA, 'N')
            my_opA = 'T';
        else
            my_opA = 'N';
        end
end
Y = mul_A_default_iter(eqn, opts, my_opA, X, 'N');
end

function Y = mfun(X, opts)
% function Y = mfun(X, opts)
% This is a function handle used specifically in pcr iterative solver
% that accepts the vector input X, and the preconditioner matrices U
% and L, to reconstruct the M preconditioner matrix, and solve the
% system Y = U \ (L \ X).

Y = opts.usfs.default_iter.PA_R \ (opts.usfs.default_iter.PA_L \ X);
end
