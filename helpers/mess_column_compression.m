function [Z, D] = mess_column_compression(Z, opZ, D, tol, info)
%	Computes a compressed representation of Z (and D).
%
%    [Z, D] = mess_column_compression(Z, opZ, D, tol, info)
%
%   Input
%       Z             matrix of interest
%
%       opZ           character specifiyng if Z should be transposed
%                       'N': Z*Z' or Z*D*Z'
%                       'T': Z'*Z or Z'*D*Z
%                     (optional, default 'N')
%
%       D             symmetric matrix of interest, if empty [] the Z*Z' or
%                     Z'*Z factorizations are considered
%                     (optional, default [])
%
%       tol           the truncation tolerance used in the rank-revealing
%                     SVD or eigenvalue decomposition
%                     (optional, default eps)
%
%       info          {0, 1}, disable/enable verbose mode
%                     (optional, default 0)
%
%   Output
%       Z             compressed low rank factor
%       D             compressed low rank factor, empty if D was empty
%                     on input

%
% This file is part of the M-M.E.S.S. project 
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2021 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%


% Check and assign input arguments.
if(issparse(Z))
    % This is just a safety measure that is hopefully never executed
    Z = full(Z);

    warning('MESS:dense',...
        ['Converting low rank factor to dense format. ' ...
        'This should never be necessary.']');
end

if nargin < 2
    opZ = 'N';
end

if strcmp(opZ, 'N')
    m = size(Z, 2);
else
    m = size(Z, 1);
end

if (nargin >= 3) && not(isempty(D))
    assert((norm(D - D', 'fro') < eps) ...
        && isequal(size(D), [m m]), ...
        'MESS:data', ...
        'The D factor has to be symmetric of size %d.', m);
else
    D = [];
end

if nargin < 4
    tol = eps;
end

if nargin < 5
    info = 0;
end

% Perform compression.
if isempty(D)
    if strcmp(opZ, 'N')
        % Z*Z' case.
        [U, S, ~] = svd(Z, 'econ');

        L = S;
        S = diag(S);
        l = length(find(S.^2 > tol * S(1)^2));
        Z = U * L(:, 1:l);

        if info
            fprintf(1, 'cc: %d -> %d  (tol: %e)\n', m, size(Z, 2), tol);
        end
    else
        % Z'*Z case.
        [~, S, V] = svd(Z, 'econ');

        L = S;
        S = diag(S);
        l = length(find(S.^2 > tol * S(1)^2));
        Z = (L(1:l, :) * V)';

        if info
            fprintf(1, 'cc: %d -> %d  (tol: %e)\n', m, size(Z, 1), tol);
        end
    end
else
    if strcmp(opZ, 'N')
        % Z*D*Z' case.
        [Q, R] = qr(Z, 0);
        RDR    = R * D * R';
        [V, S] = eig(0.5 * (RDR + RDR'));
        S      = diag(S);

        r = abs(S) > tol * max(abs(S));
        Z = Q * V(:, r);
        D = diag(S(r));

        if info
            fprintf(1, 'cc: %d -> %d  (tol: %e)\n', m, size(Z, 2), tol);
        end
    else
        % Z'*D*Z case.
        [Q, R] = qr(Z', 0);
        RDR    = R * D * R';
        [V, S] = eig(0.5 * (RDR + RDR'));
        S      = diag(S);

        r = abs(S) > tol * max(abs(S));
        Z = (Q * V(:, r))';
        D = diag(S(r));

        if info
            fprintf(1, 'cc: %d -> %d  (tol: %e)\n', m, size(Z, 1), tol);
        end
    end
end
