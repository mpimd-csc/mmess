function X = sol_E_so_iter(eqn, opts, opE, B, opB)

% function X=sol_E_so_iter(eqn, opts, opE, B, opB)
%
% This function returns X = E_\B, where matrix E_ given by
% structure eqn and input matrix B could be transposed. Matrix E_
% is assumed to be quadratic and has the same size as A_ in
% structure eqn.
%
%   Inputs:
%
%   eqn       structure containing field 'E_'
%   opts      structure containing parameters for the algorithm.
%   opE       character specifying the shape of E_
%                  opE = 'N' solves E_*X = opB(B)
%                  opE = 'T' solves E_'*X = opB(B)
%   B         p-x-q matrix
%   opB       character specifying the shape of B
%                  opB = 'N' solves  opE(E_)*X = B
%                  opB = 'T' solves  opE(E_)*X = B'
%
%   Output:
%
%   X         matrix fulfilling equation  opE(E_)*X = opB(B)
%
% This function uses another function size_so_iter(eqn,
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
if not(ischar(opE)) || not(ischar(opB))
    mess_err(opts, 'error_arguments', 'opE or opB is not a char');
end
opE = upper(opE);
opB = upper(opB);
if not(opE == 'N' || opE == 'T')
    mess_err(opts, 'error_arguments', 'opE is not ''N'' or ''T''');
end

if not(opB == 'N' || opB == 'T')
    mess_err(opts, 'error_arguments', 'opB is not ''N'' or ''T''');
end

if (not(isnumeric(B))) || (not(ismatrix(B)))
    mess_err(opts, 'error_arguments', 'B has to be a matrix');
end

% Initial guess for vector X

if not(isfield(opts.usfs.so_iter, 'X0_E'))
    opts.usfs.so_iter.X0_E = [];
end

%% Check data in eqn structure

switch lower(opts.usfs.so_iter.method_E)
    case 'pcg'
        if not(opts.usfs.so_iter.E_is_spd)
            mess_err(opts, 'error_arguments', ...
                     'Field eqn.E_ is not symmetric and positive definite');
        end
    case {'minres', 'symmlq', 'pcr'}
        if not(issymmetric(eqn.E_))
            mess_err(opts, 'error_arguments', 'Field eqn.E_ is not symmetric');
        end
end

n = 2 * size_so_iter(eqn, opts);

%% Preallocate solution
if opB == 'N'
    X = zeros(size(B));
    flags = zeros(1, size(B, 2));
else
    X = zeros(size(B'));
    flags = zeros(1, size(B, 1));
end

%% Create anonymous functions
% To call multiplication with E respecting opE

switch lower(opts.usfs.so_iter.method_E)

    case {'bicg', 'lsqr', 'qmr'}
        mul_E = @(X, flag) flagged_mul_E(X, flag, eqn, opts, opE);

    otherwise
        mul_E = @(X) mul_E_so_iter(eqn, opts, opE, X, 'N');

end

% For calling the actual iterative solver
solver = eval(sprintf('@%s', lower(opts.usfs.so_iter.method_E)));

%% Perform solve operations
switch opB

    case 'N'
        if not(n == size(B, 1))
            mess_err(opts, 'error_arguments', ...
                     ['Number of rows of E_ differs from number ' ...
                      'of rows of B']);
        end

        switch lower(opts.usfs.so_iter.method_E)

            case 'gmres'
                for i = 1:size(B, 2)
                    [X(:, i), flags(i)] = ...
                        solver(mul_E, B(:, i), ...
                               opts.usfs.so_iter.restIter, ...
                               opts.usfs.so_iter.res_tol, ...
                               opts.usfs.so_iter.max_iter, ...
                               opts.usfs.so_iter.PE_L, ...
                               opts.usfs.so_iter.PE_R, ...
                               opts.usfs.so_iter.X0_E);
                end
            case 'pcr'
                for i = 1:size(B, 2)
                    [X(:, i), flags(i)] = ...
                        solver(mul_E, B(:, i), ...
                               opts.usfs.so_iter.res_tol, ...
                               opts.usfs.so_iter.max_iter, ...
                               @(X) mfun(X, opts), ...
                               opts.usfs.so_iter.X0_E);
                end
            otherwise
                for i = 1:size(B, 2)
                    [X(:, i), flags(i)] = ...
                        solver(mul_E, B(:, i), ...
                               opts.usfs.so_iter.res_tol, ...
                               opts.usfs.so_iter.max_iter, ...
                               opts.usfs.so_iter.PE_L, ...
                               opts.usfs.so_iter.PE_R, ...
                               opts.usfs.so_iter.X0_E);
                end

        end

    case 'T'
        if not(n == size(B, 2))
            mess_err(opts, 'error_arguments', ...
                     ['Number of rows of E_ differs from number ' ...
                      'of columns of B']);
        end

        switch lower(opts.usfs.so_iter.method_E)

            case 'gmres'
                for i = 1:size(B, 1)
                    [X(:, i), flags(i)] = ...
                        solver(mul_E, B(i, :)', ...
                               opts.usfs.so_iter.restIter, ...
                               opts.usfs.so_iter.res_tol, ...
                               opts.usfs.so_iter.max_iter, ...
                               opts.usfs.so_iter.PE_L, ...
                               opts.usfs.so_iter.PE_R, ...
                               opts.usfs.so_iter.X0_E);
                end
            case 'pcr'
                for i = 1:size(B, 1)
                    [X(:, i), flags(i)] = ...
                        solver(mul_E, B(i, :)', ...
                               opts.usfs.so_iter.res_tol, ...
                               opts.usfs.so_iter.max_iter, ...
                               @(X) mfun(X, opts), ...
                               opts.usfs.so_iter.X0_E);
                end
            otherwise
                for i = 1:size(B, 1)
                    [X(:, i), flags(i)] = ...
                        solver(mul_E, B(i, :)', ...
                               opts.usfs.so_iter.res_tol, ...
                               opts.usfs.so_iter.max_iter, ...
                               opts.usfs.so_iter.PE_L, ...
                               opts.usfs.so_iter.PE_R, ...
                               opts.usfs.so_iter.X0_E);
                end

        end
end
if any(flags)
    mess_warn(opts, 'usfs_iter', ...
              [lower(opts.usfs.so_iter.method_E) ...
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

function Y = flagged_mul_E(X, flag, eqn, opts, opE)
% function Y = flagged_mul_E(X, flag, eqn, opts, opE)
% This is a function handle that accepts the vector input X and the matrix
% E_ given by the eqn structure, and returns the matrix vector product
% A_*X. The input 'flag' defines whether E_ should be transposed or not.
switch lower(flag)
    case 'notransp'
        my_opE = opE;
    case 'transp'
        if strcmp(opE, 'N')
            my_opE = 'T';
        else
            my_opE = 'N';
        end
end
Y = mul_E_so_iter(eqn, opts, my_opE, X, 'N');
end

function Y = mfun(X, opts)
% function Y = mfun(X, opts)
% This is a function handle used specifically in pcr iterative solver
% that accepts the vector input X, and the preconditioner matrices U
% and L, to reconstruct the M preconditioner matrix, and solve the
% system Y = U \ (L \ X).

Y = opts.usfs.so_iter.PE_R \ (opts.usfs.so_iter.PE_L \ X);
end
