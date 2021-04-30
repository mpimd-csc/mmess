function weights = compute_quadrature_weights(h, nodes)
% Compute the weights for a quadrature formula on [0,h] defined by the
% given nodes in the interval. If length(h) == 2 then use the interval
% [h(1) h(2)] instead.

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2021 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%


    Cf = fliplr(vander(nodes))';
    if length(h) > 1
        b = h(2).^(1:length(nodes))' ...
            ./ (1:length(nodes))' - h(1).^(1:length(nodes))' ...
            ./ (1:length(nodes))';
    else
        b = h.^(1:length(nodes))' ./ (1:length(nodes))';
    end
%     weights = Cf\b;
    % Better for high number of nodes (ill-conditioned system)
    weights = vandersolve(Cf, b);


    function x = vandersolve(M,b)
    % Solve the Vandermonde system Mx = b according to the method described
    % in \r{A}ke Bj\"{o}rck and Victor Pereyra, "Solution of Vandermonde
    % Systems of Equations", Mathematics of Computation, Vol. 24, No. 112
    % (Oct., 1970), pp. 893-903, http://www.jstor.org/stable/2004623.
    % (Eq. (14)-(15) for the primal system.)

        n = size(M,1) - 1;
        d = b;
        for k=1:n
            d(k+1:n+1) = d(k+1:n+1) - M(2,k)*d(k:n);
        end

        x = d;
        for k=n:-1:1
            x(k+1:n+1) = x(k+1:n+1) ./ (M(2,k+1:n+1) - M(2,1:n-k+1))';
            x(k:n) = x(k:n) - x(k+1:n+1);
        end

    end

end