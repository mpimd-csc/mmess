function X = sol_ApE_default_iter(eqn, opts, opA, p, opE, C, opC)

% function X=sol_ApE_default_iter(eqn, opts,opA,p,opE,C,opC)
%
% This function returns X = (A_ + p*E_)\C, where matrices A_ and E_
% given by structure eqn and input matrix C could be transposed.
% Matrices A_ and E_ are assumed to be quadratic.
%
%   Inputs:
%
%   eqn     structure containing fields 'A_' and 'E_'
%   opts    structure containing parameters for the algorithm
%   opA     character specifying the shape of A
%           opA = 'N' solves (A_ + p* opE(E_))*X = opC(C)
%           opA = 'T' solves (A_' + p* opE(E_))*X = opC(C)
%   p       scalar value
%   opE     character specifying the shape of E_
%           opE = 'N' solves (opA(A_) + p* E_)*X = opC(C)
%           opE = 'T' solves (opA(A_) + p* E_')*X = opC(C)
%   C       n-x-p matrix
%   opC     character specifies the form of opC(C)
%           opC = 'N' solves (opA(A_) + p* opE(E_))*X = C
%           opC = 'T' solves (opA(A_) + p* opE(E_))*X = C'
%
%   Output:
%
%   X       matrix fulfilling equation (opA(A_)+p*opE(E_))*X = opC(C)
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

if (not(isnumeric(p))) || not(length(p) == 1)
    mess_err(opts, 'error_arguments', 'p is not numeric');
end

if (not(isnumeric(C))) || (not(ismatrix(C)))
    mess_err(opts, 'error_arguments', 'C has to be a matrix');
end

n = size_default_iter(eqn, opts);

% Initial guess for vector X

if not(isfield(opts.usfs.default_iter, 'X0_ApE'))
    opts.usfs.default_iter.X0_ApE = [];
end
%% Check data in eqn structure
switch lower(opts.usfs.default_iter.method_ApE)
    case {'minres', 'pcg', 'symmlq'}
        if not(issymmetric(eqn.A_ + p * eqn.E_))
            mess_err(opts, 'error_arguments', 'Resulting matrix of (eqn.A_+p*eqn.E_) is not symmetric');
        end
end
%% Preallocate solution

if opC == 'N'
    X = zeros(size(C));
    flags = zeros(1, size(C, 2));
else
    X = zeros(size(C'));
    flags = zeros(1, size(C, 1));
end

%% Create anonymous functions
% To call multiplication with ApE respecting opA and opE
switch lower(opts.usfs.default_iter.method_ApE)

    case {'bicg', 'lsqr', 'qmr'}

        mul_ApE = @(X, flag) flagged_mul_ApE(X, flag, eqn, opts, opA, p, opE);
    case 'pcg'

        mul_ApE = @(X) -mul_ApE_default_iter(eqn, opts, opA, p, opE, X, 'N');
    otherwise

        mul_ApE = @(X) mul_ApE_default_iter(eqn, opts, opA, p, opE, X, 'N');
end

% For calling the actual iterative solver
solver = eval(sprintf('@%s', lower(opts.usfs.default_iter.method_ApE)));
%% Perform solve operations

switch opC

    case 'N'
        if not(n == size(C, 1))
            mess_err(opts, 'error_arguments', ['Number of rows ' ...
                                               'of A_ differs with number ' ...
                                               'of rows of C']);
        end

        switch lower(opts.usfs.default_iter.method_ApE)

            case 'pcg'
                for i = 1:size(C, 2)
                    [x, flags(i)] = ...
                        solver(mul_ApE, C(:, i), ...
                               opts.usfs.default_iter.res_tol, ...
                               opts.usfs.default_iter.max_iter, ...
                               opts.usfs.default_iter.PApE_L, ...
                               opts.usfs.default_iter.PApE_R, ...
                               opts.usfs.default_iter.X0_ApE);
                    X(:, i) = -x;
                end

            case 'gmres'
                for i = 1:size(C, 2)
                    [X(:, i), flags(i)] = ...
                        solver(mul_ApE, C(:, i), ...
                               opts.usfs.default_iter.restIter, ...
                               opts.usfs.default_iter.res_tol, ...
                               opts.usfs.default_iter.max_iter, ...
                               opts.usfs.default_iter.PApE_L, ...
                               opts.usfs.default_iter.PApE_R, ...
                               opts.usfs.default_iter.X0_ApE);
                end
            case 'pcr'
                for i = 1:size(C, 2)
                    [X(:, i), flags(i)] = ...
                        solver(mul_ApE, C(:, i), ...
                               opts.usfs.default_iter.res_tol, ...
                               opts.usfs.default_iter.max_iter, ...
                               @(X) mfun(X, opts), ...
                               opts.usfs.default_iter.X0_ApE);
                end
            otherwise
                for i = 1:size(C, 2)
                    [X(:, i), flags(i)] = ...
                        solver(mul_ApE, C(:, i), ...
                               opts.usfs.default_iter.res_tol, ...
                               opts.usfs.default_iter.max_iter, ...
                               opts.usfs.default_iter.PApE_L, ...
                               opts.usfs.default_iter.PApE_R, ...
                               opts.usfs.default_iter.X0_ApE);
                end

        end

    case 'T'
        if not(n == size(C, 2))
            mess_err(opts, 'error_arguments', ...
                     ['Number of rows of A_ differs with number ' ...
                      'of columns of C']);
        end

        switch lower(opts.usfs.default_iter.method_ApE)

            case 'pcg'
                for i = 1:size(C, 1)
                    [x, flags(i)] = ...
                        solver(mul_ApE, C(i, :)', ...
                               opts.usfs.default_iter.res_tol, ...
                               opts.usfs.default_iter.max_iter, ...
                               opts.usfs.default_iter.PApE_L, ...
                               opts.usfs.default_iter.PApE_R, ...
                               opts.usfs.default_iter.X0_ApE);
                    X(:, i) = -x;
                end

            case 'gmres'
                for i = 1:size(C, 1)
                    [X(:, i), flags(i)] = ...
                        solver(mul_ApE, C(i, :)', ...
                               opts.usfs.default_iter.restIter, ...
                               opts.usfs.default_iter.res_tol, ...
                               opts.usfs.default_iter.max_iter, ...
                               opts.usfs.default_iter.PApE_L, ...
                               opts.usfs.default_iter.PApE_R, ...
                               opts.usfs.default_iter.X0_ApE);
                end
            case 'pcr'
                for i = 1:size(C, 1)
                    [X(:, i), flags(i)] = ...
                        solver(mul_ApE, C(i, :)', ...
                               opts.usfs.default_iter.res_tol, ...
                               opts.usfs.default_iter.max_iter, ...
                               @(X) mfun(X, opts), ...
                               opts.usfs.default_iter.X0_ApE);
                end
            otherwise
                for i = 1:size(C, 1)
                    [X(:, i), flags(i)] = ...
                        solver(mul_ApE, C(i, :)', ...
                               opts.usfs.default_iter.res_tol, ...
                               opts.usfs.default_iter.max_iter, ...
                               opts.usfs.default_iter.PApE_L, ...
                               opts.usfs.default_iter.PApE_R, ...
                               opts.usfs.default_iter.X0_ApE);
                end
        end

end
if any(flags)
    mess_warn(opts, 'usfs_iter', ...
              [lower(opts.usfs.default_iter.method_ApE) ...
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

function Y = flagged_mul_ApE(X, flag, eqn, opts, opA, p, opE)
% function Y = flagged_mul_ApE(X, flag, eqn, opts, opA, p, opE)
% This is a function handle that accepts the vector input X, the scalar p
% and the matrices A_ and E_ given by the eqn structure, and returns the matrix
% vector product (A_+p*E_)*X. The input 'flag' defines whether A_ and E_
% should be transposed or not.

switch lower(flag)
    case 'notransp'
        Y = mul_ApE_default_iter(eqn, opts, opA, p, opE, X, 'N');
    case 'transp'
        if strcmp(opA, 'N')
            my_opA = 'T';
        else
            my_opA = 'N';
        end
        if strcmp(opE, 'N')
            my_opE = 'T';
        else
            my_opE = 'N';
        end
        Y = mul_ApE_default_iter( ...
                                 eqn, opts, my_opA, conj(p), my_opE, X, 'N');
end

end

function Y = mfun(X, opts)
% function Y = mfun(X, opts)
% This is a function handle used specifically in pcr iterative solver
% that accepts the vector input X, and the preconditioner matrices U
% and L, to reconstruct the M preconditioner matrix, and solve the
% system Y = U \ (L \ X).

Y = opts.usfs.default_iter.PApE_R \ (opts.usfs.default_iter.PApE_L \ X);
end
