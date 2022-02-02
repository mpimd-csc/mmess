function [gamma] = compute_additive_coefficients(order, symmetric)
% Compute order conditions for the asymmetric and symmetric additive
% splitting schemes
%
% The order must be even for the symmetric schemes.
%
% Use symmetric = true for the symmetric version, symmetric = false for the
% asymmetric.
%

%
% This file is part of the M-M.E.S.S. project 
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

if symmetric
    s = order / 2;

    A = zeros(s,s);
    A(1,:) = 1;
    for k = 1:s-1
        A(k+1, :) = (1:s).^(-2*k);
    end

    b = [1/2; zeros(s-1, 1)];

    gamma = A\b;
else
    s = order;

    A = zeros(s,s);
    A(1,:) = 1;
    for k = 1:s-1
        A(k+1, :) = (1:s).^(-k);
    end

    b = [1; zeros(s-1, 1)];

    gamma = A\b;
end