function [eqn, opts, oper] = KSM_compute_shifts(eqn, opts, oper)
% function [eqn,opts,oper] = KSM_compute_shifts(eqn,opts,oper)
% Function that computes the new shift of the rational Krylov subspace
%
% Input and output:
%
%   eqn       structure containing equation data
%
%   opts      structure containing parameters for the algorithm
%
%   oper      contains function handles with operations for A and E
%
%
% Output:
%
%  opts.KSM.compute_struct       structure containing the following:
%
%    new_shift    computed shift to use in the basis construction of RK
%
%    shifts   all the computed shifts
%
%    shift_cmplxflag   true if the last computed shift is complex.
%                      Then, we need to incorporate its complex conjugate
%

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

p = size(eqn.W, 2);

it = opts.KSM.compute_struct.it;
if it == 1

    % compute the eigenvalues of T
    eH = eig(opts.KSM.compute_struct.T(1:it * p, 1:it * p));
    eHpoints = sort(opts.KSM.init_shifts);

    % compute the new shift
    opts.KSM.compute_struct.new_shift = ...
        newpolei(eHpoints, eH, opts.KSM.init_shifts(2) * ones(p, 1));
    % check that it's contained in the right half plane
    if real(opts.KSM.compute_struct.new_shift) < 0
        opts.KSM.compute_struct.new_shift = ...
            -real(opts.KSM.compute_struct.new_shift) + ...
            1i * imag(opts.KSM.compute_struct.new_shift);
    end
    % update the shift array
    opts.KSM.compute_struct.shifts = ...
        [opts.KSM.compute_struct.shifts, ...
         opts.KSM.compute_struct.new_shift];

    % If pole is complex, include its conjugate
    opts.KSM.compute_struct.shift_cmplxflag = ...
        not(isreal(opts.KSM.compute_struct.new_shift));

else

    if strcmp(opts.KSM.type_eqn, 'CARE') && ...
            isfield(opts.KSM.compute_struct, 'Y')
        if isfield(opts.KSM, 'CARE_shifts') && ...
                strcmp(opts.KSM.CARE_shifts, 'Ritz_closedloop')
            Y = opts.KSM.compute_struct.Y;
            B = opts.KSM.compute_struct.Bm;
            B = B(1:p * it, :);
            sizeY = size(Y);
            % since we may not compute Y at each iteration, we pad it with
            % extra zero rows and columns to match the dimension of T if
            % necessary
            if not(sizeY(1) == p * it)
                Y = [[Y; zeros(it * p - sizeY(1), sizeY(2))], ...
                     zeros(it * p, it * p - sizeY(2))];
            end
            eH = ...
                eig(opts.KSM.compute_struct.T(1:it * p, 1:it * p) + ...
                    B * (B' * Y));
        else
            eH = eig(opts.KSM.compute_struct.T(1:it * p, 1:it * p));
        end
    else
        eH = eig(opts.KSM.compute_struct.T(1:it * p, 1:it * p));
    end
    eH = sort(eH);
    eHorig = eH;

    if strcmp(opts.KSM.type_shifts, 'complex')
        % Complex poles. Compute set for next complex pole of r_m

        if any(imag(eH)) && max(abs(imag(eH))) > 1e-5 && length(eH) > 2
            % Roots lambdas come from convex hull too
            eH = [eH; -opts.KSM.init_shifts(1)];
            ij = convhull(real(eH), imag(eH));
            eH = eH(ij);
            ieH = length(eH);
            missing = it * p - ieH;
            while missing > 0
                % include enough points from the border
                neweH = (eH(1:ieH - 1) + eH(2:ieH)) / 2;
                eH = [eH; neweH]; %#ok<AGROW>
                missing = it * p - length(eH);
            end

            eHpoints = -eH;
            eH = eHorig;
        else
            % if all real eigs, no convex hull possible
            eHpoints = sort( ...
                            [opts.KSM.init_shifts(2); ...
                             opts.KSM.init_shifts(1).'; ...
                             -real(eH)]);
        end

    else
        % Real poles s from real set. Compute complex roots of r_m via Ritz
        % convex hull
        if any(imag(eH)) && length(eH) > 2
            % Roots lambdas come from convex hull too
            eH = [eH; -opts.KSM.init_shifts(2); -opts.KSM.init_shifts(1).'];
            ij = convhull(real(eH), imag(eH));
            eH = eH(ij);
            ieH = length(eH);
            missing = it * p - ieH;
            while missing > 0
                % include enough points from the border
                neweH = (eH(1:ieH - 1) + eH(2:ieH)) / 2;
                eH = [eH; neweH]; %#ok<AGROW>
                missing = it * p - length(eH);
            end
            eH = eH(1:it * p);
        end
        eHpoints = sort( ...
                        [opts.KSM.init_shifts(2); ...
                         opts.KSM.init_shifts(1).'; ...
                         -real(eH)]);
        eH = eHorig;
    end

    gs = kron(opts.KSM.compute_struct.shifts(2:end), ones(1, p))';
    opts.KSM.compute_struct.new_shift = newpolei(eHpoints, eH, gs);
    % check that it's contained in the right half plane
    if real(opts.KSM.compute_struct.new_shift) < 0
        opts.KSM.compute_struct.new_shift = ...
            -real(opts.KSM.compute_struct.new_shift) + ...
            1i * imag(opts.KSM.compute_struct.new_shift);
    end

    % If pole is complex, include its conjugate
    opts.KSM.compute_struct.shift_cmplxflag = ...
        not(isreal(opts.KSM.compute_struct.new_shift));

    % update the shift array
    opts.KSM.compute_struct.shifts = ...
        [opts.KSM.compute_struct.shifts, ...
         opts.KSM.compute_struct.new_shift];

end
% print the shift if needed
if opts.shifts.info
    mess_fprintf(opts, 'New shift: %10.5e\n', ...
                 opts.KSM.compute_struct.new_shift);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r = ratfun(x, eH, s)
r = zeros(length(x), 1);
for j = 1:length(x)
    r(j) = abs(prod((x(j) - s) ./ (x(j) - eH)));
end

return

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function snew = newpolei(eHpoints, eH, s)

snew_vec = zeros(length(eHpoints) - 1, 1);
for j = 1:length(eHpoints) - 1
    sval = linspace(eHpoints(j), eHpoints(j + 1), 20);
    [~, jx] = max (abs(ratfun(sval, eH, s)));
    snew_vec(j) = sval(jx);
end
[~, jx] = max(abs(ratfun(snew_vec, eH, s)));
snew = snew_vec(jx);
return
