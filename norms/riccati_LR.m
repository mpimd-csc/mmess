function [res] = riccati_LR(W, DeltaK, opts, S, S_K)
% RICCATI_LR Norm of the Riccati residual
%   Using low-rank formulation
%   R(X) = U * D * U^T
%   Since U^T * U * D is not symmetric ||U * D * U^T|| is not the
%   same as ||U^T * U * D||, but the non-zero eigenvalues of (U * D
%   * U^T) and (U^T * U * D) are equal.

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%% Check input data

if not(isfield(opts, 'bdf'))
    opts.bdf = [];
end

if isstruct(opts.bdf) && isfield(opts.bdf, 'tau') && ...
        isfield(opts.bdf, 'beta') && isempty(S_K)

    DeltaK = DeltaK * sqrt(opts.bdf.tau * opts.bdf.beta);
    % DeltaK needs to be scaled with tb to compute the residual norm of the
    % ARE; DeltaK does not contain tb since it's for the DRE; if S_K is not
    % empty DeltaK results from a line search and also contains parts from
    % W; the scaling with tb is done in mess_lrnm already.
end

U = [W, DeltaK];

if opts.LDL_T
    if isempty(S_K)
        D = blkdiag(S, -eye(size(DeltaK, 2)));
    else
        D = blkdiag(S, -S_K);
    end
else
    D = blkdiag(eye(size(W, 2)), ...
                -eye(size(DeltaK, 2)));
end

%% Compute norm
if opts.norm == 2 % 2-norm
    f = @(x) max(abs(eig(x)));
elseif strcmpi(opts.norm, 'fro') % Frobenius norm
    f = @(x) norm(eig(x));
else
    mess_err(opts, 'riccati_LR: Unsupported norm');
end

res = f(U' * U * D);

end
