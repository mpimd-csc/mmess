function [ opts, out, stop ] = prepare_next_adi_iteration( opts, out, res, rc, outer_res, i)
% Evaluate stopping criteria of LRADI for exact and inexact case and check
% whether line search is necessary.

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%


%% Check input

%% Evaluate stopping criteria
stop = 0;

% Riccati tolerance reached, stop
if not(isempty(outer_res)) && ...
   (outer_res(i) < opts.nm.res_tol)

    stop = 1;

% Inexact ADI
elseif opts.adi.inexact

    % Lyapunov residual is growing, stop and perform line search
    if opts.adi.res_tol && ...
       (i > 2) && ...
       (res(i) > res(1))

        out.linesearch = 1;
        stop = 1;

    % Outer tolerance reached and not close to finish Newton iteration
    elseif opts.adi.res_tol && ...
           (res(i) < opts.adi.outer_tol) && ...
           ((res(i) >= 1e2 * opts.nm.res_tol) || ...
            ((i > 1) && ...
            (outer_res(i) >= outer_res(i - 1))))
        
        stop = 1;

    % else
    % Outer tolerance not reached and Lyapunov residual is not growing or
    % outer tolerance is reached but Newton iteration is almost finished,
    % do NOT stop ADI iteration!
    end

% ADI tolerance reached, stop ADI iteration
elseif ((opts.adi.res_tol && ...
        (res(i) < opts.adi.res_tol)) || ...
        (opts.adi.rel_diff_tol && ...
         (rc(i) < opts.adi.rel_diff_tol))) && ...
       (isempty(outer_res) || ... 
        ((res(i) >= 1e2 * opts.nm.res_tol) || ...
         ((i > 1) && ...
          (outer_res(i) >= outer_res(i - 1)))))

    stop = 1;

    % Riccati residual is growing, perform line search
    if not(isempty(outer_res)) && ...
       outer_res(i) > outer_res(1)
        
        out.linesearch = 1;
    end

% Lyapunov residual is growing and inexact ADI with line search probably
% failed already, stop ADI iteration and restart Newton iteration with
% exact ADI iteration
elseif opts.adi.res_tol && ...
       (res(i) > res(1) * 1e2) && ...
       isfield(opts,'nm')

    out.restart = 1;
    stop = 1;
end
end

