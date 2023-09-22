function [eqn, opts, oper, nShifts] = ...
    mess_lrradi_get_shifts_hamOpti_generalized(eqn, opts, oper, W, Z, Y)
% Compute the Hamiltonian residual shifts for the RADI method.
% Optionally, run residual minimizing search afterwards, and return
% optimized shifts.
%
% The routine avoids solving linear systems with E. The shifts are computed
% as generalized Ritz values of the matrix pair (H, EE), where
%    H = [A + U*V'  B*B'         ], EE = [ E 0  ]
%         W'*W      -(A + U*V')' ]       [ 0 E' ].
%
% The matrix pair is projected onto the subspace spanned by the ell last vectors
% in Z. The generalized eigenvalues and eigenvectors of the obtained pair
% (Hproj, Eproj) are then used to compute the shift in the following way:
%     - all obtained eigenvectors [x; y] are sorted by the value of
%       tau = norm(y)^2 / abs( x'*E'*y ).
%     - the eigenvalue corresponding to the largest value of tau is chosen
%       as the next shift.
%
% Input
%   eqn                 struct contains data for equations
%
%   opts                struct contains parameters for the algorithm
%
%   oper                struct contains function handles for operation
%                       with A and E
%
%   W                   the current residual matrix in RADI
%
%   Z                   the current Z matrix in RADI
%
%   Y                   the small square factor in the RADI solution
%                       (not used)
%
%
% Output
%   opts                the modified input structure
%
%   nShifts             the number of shifts computed (=1, always)
%
%
% Input fields in struct eqn used by this routine:
%   eqn.BB      dense (n x m1) matrix from quadratic term
%
%   eqn.U       dense (n x m3) matrix U from rank-m3 Update
%
%   eqn.V       dense (n x m3) matrix V from rank-m3 Update
%
%   eqn.type    possible  values: 'N', 'T'
%               determining whether (N) or (T) is solved
%
%   eqn.haveE   possible  values: false, true
%               if haveE = false: matrix E in eqn.E_ is assumed to be identity
%
%
% Input fields in struct opts used by this routine:
%   shifts.history
%       Possible values: inf or k*p, where k>=2, and p = size(W, 2).
%       The number of the last columns of Z which are used to project the
%       Hamiltonian pair, i.e., the dimension of the projection subspace.
%       If set to inf, all columns of Z will be used.
%
%   shifts.naive_update_mode
%       Possible values: true | [false]
%       Sets the method of updating the projection between successive calls
%       to mess_lrradi_get_shifts_hamOpti_generalized. If set to true, the QR
%       factorization will be recomputed from scratch in each call to
%       obtain the orthonormal basis for the projection subspace. If set
%       to (default) false, the QR factorization from the previous step
%       will be updated to save computation time.
%
%
% Output fields in struct opts set by this routine:
%   shifts.p
%       Contains the computed shift.
%
%   opts.shifts.tmp
%       struct with all the temporary matrices for the shifts
%

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

% Determine the number of inputs.
p = size(W, 2);

% Check if the projection subspace size is valid.
if not(opts.shifts.history == Inf || ...
       (mod(opts.shifts.history, p) == 0) && ...
       (opts.shifts.history >= 2 * p))
    mess_warn(opts, 'control_data', ['Shifts history length should be ', ...
                                     'Inf, or a multiple of number of ', ...
                                     'residual columns p (>=2*p).', ...
                                     'Set opts.shifts.history = %i.'], ...
              6 * p);
    opts.shifts.history = 6 * p;
end

if not(isfield(opts.shifts, 'naive_update_mode')) || ...
        not(islogical(opts.shifts.naive_update_mode))
    opts.shifts.naive_update_mode = false;
end

% Initialize the shift generator when this function is called for the
% first time.
if not(isfield(opts.shifts, 'tmp')) || not(isfield(opts.shifts.tmp, 'U'))
    [eqn, opts, oper, W, Z, Y] = initialize(eqn, opts, oper, W, Z, Y);
end

% Compute the shift. Truncate shift to real if the imaginary part is small
% enough.
[eqn, opts, oper, shift] = getNextShift(eqn, opts, oper, W, Z, Y);
if (abs(imag(shift)) / abs(shift)) < 1e-8
    shift = real(shift);
end

opts.shifts.p = shift;
nShifts       = 1;

end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliary functions for the shift generator.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [eqn, opts, oper, W, Z, Y] = initialize(eqn, opts, oper, W, Z, Y)
% Queue for storing previous subspace dimensions.
opts.shifts.tmp.block_queue        = {};
opts.shifts.tmp.dim_current        = 0;
opts.shifts.tmp.dim_total_subspace = 0;
opts.shifts.tmp.dim_max            = opts.shifts.history; % .. Alias ..

if isfield(eqn, 'Z0') && (size(eqn.Z0, 2) > opts.shifts.history)
    % Ensure that the shift history is larger than the initial
    % solution.
    opts.shifts.tmp.dim_max = opts.shifts.tmp.dim_max + size(eqn.Z0, 2);
end
end

function [eqn, opts, oper, shift] = getNextShift(eqn, opts, oper, W, Z, Y)
% Computes the next shift.

% Update (or initialize) orthogonal basis of the projected subspace.
[eqn, opts, oper, W, Z, Y] = update_qr(eqn, opts, oper, W, Z, Y);

% Assemble/update the projected Hamiltonian matrix (size = 2*ell).
[eqn, opts, ~, ~, ~, ~] = update_ray(eqn, opts, oper, W, Z, Y);

% Compute the eigenpairs of the pair (Hproj, Eproj).
if eqn.haveE
    % Compute the eigenvalues of the matrix pencil (Hproj, Eproj).
    % It turns out that reducing to ordinary eigenvalue problem is
    % faster.
    % [XX, LL] = eig(opts.shifts.tmp.Hproj, opts.shifts.tmp.Eproj);
    [XX, LL] = eig(opts.shifts.tmp.Eproj \ opts.shifts.tmp.Hproj);
else
    [XX, LL] = eig(opts.shifts.tmp.Hproj);
end

% Filter out the stable eigenvalues.
LL   = diag(LL);
perm = (real(LL) < 0);
LL   = LL(perm);
XX   = XX(:, perm);

% Rescale the stable eigenvectors.
nX = length(XX) / 2;

% Find the eigenvector with the largest norm(y)^2 / abs(xt*Et*y)
taj = 1;
maxi = -1;

for i = 1:length(LL)
    temp = XX(nX + 1:end, i); % y
    curr = norm(temp);      % norm(y)

    if eqn.haveE
        temp = opts.shifts.tmp.UtEU' * temp; % E' * y
    end

    temp = abs(XX(1:nX, i)' * temp); % x'*E'*y
    if curr
        curr = curr / temp * curr; % norm(y)^2 / abs(x'*E'*y)
    end

    if (i == 1) || (curr > maxi)
        maxi = curr;
        taj  = i;
    end
end

% The shift is the eigenvalue associated to taj.
shift = LL(taj);
end

function [eqn, opts, oper, W, Z, Y] = update_qr(eqn, opts, oper, W, Z, Y)
% Updates the orthogonal basis for the projection subspace.
dim_current        = opts.shifts.tmp.dim_current;
dim_total_subspace = opts.shifts.tmp.dim_total_subspace;
dim_max            = opts.shifts.tmp.dim_max;

% Detect how many columns have been added, and how many have to be
% dropped from the previous step.
k       = size(Z, 2);
dim_new = k - dim_total_subspace;
opts.shifts.tmp.block_queue{length(opts.shifts.tmp.block_queue) + 1} = dim_new;

dim_drop = 0;
dim_keep = dim_current;
while dim_keep + dim_new > dim_max
    dim_block                      = opts.shifts.tmp.block_queue{1};
    opts.shifts.tmp.block_queue(1) = [];

    % Drop columns block by block from the queue.
    dim_keep = dim_keep - dim_block;
    dim_drop = dim_drop + dim_block;
end

% If the first shift is required, use the orthonormal basis of W.
% This is NOT reused for the following shifts.
if k == 0
    [opts.shifts.tmp.U, opts.shifts.tmp.RR] = qr(W, 0);

    return
end

% Otherwise, update the QR factorization.
if opts.shifts.naive_update_mode || not(isfield(opts.shifts.tmp, 'RR'))
    % Just compute the QR factorization of the last
    % dim_keep + dim_new columns of Z.
    % May be inefficient!
    [opts.shifts.tmp.U, opts.shifts.tmp.RR] = ...
        qr(Z(:, k - (dim_keep + dim_new) + 1:k), 0);
else
    % Update the QR factorization by dropping first dim_drop columns
    % from the factorization.

    % First, compute the QRF of the last dim_keep columns of RR
    % (This could actually be done by a sequence of small Householder
    % reflectors, but it is more efficient?)
    [opts.shifts.tmp.Q, opts.shifts.tmp.RR] = ...
        qr(opts.shifts.tmp.RR(:, ...
                              dim_current - dim_keep + 1:dim_current), 0);

    % Q = Ortho. matrix that transforms old U ([drop keep]) into new
    % U ([keep]). Update the basis U with Q from the right.
    opts.shifts.tmp.U = opts.shifts.tmp.U * opts.shifts.tmp.Q;

    % The new subspace to append the basis.
    U_new = Z(:, k - dim_new + 1:k);

    % Check if a part of the new subspace is already a part of the
    % old one.
    if (size(opts.shifts.tmp.U, 2) > 0) && (size(U_new, 2) > 0)
        [U_new, ~] = qr(U_new, 0);
        [~, s, v]  = svd(opts.shifts.tmp.U' * U_new);
        if (size(s, 2) > 1) && (size(s, 1) > 1)
            % Tolerance on the cosine of the canonical angles,
            % hardcoded...
            perm = (diag(s) < 1 - 1e-11);
        else
            % A weird special case. Matlab would make diag(s) a matrix
            % here.
            perm = (s(1) < 1 - 1e-11);
        end
        % Keep only the part of the subspace already not contained
        % in U.
        U_new   = U_new * v(:, perm);
        dim_new = size(U_new, 2);
    end

    opts.shifts.tmp.RR = blkdiag(opts.shifts.tmp.RR, zeros(dim_new));

    % Orthogonalize columns of U_new against U, twice is enough.
    p = size(opts.shifts.tmp.U, 2); % do it all at once
    jj_start = 1;
    jj_end   = dim_new;
    for twice = 1:2
        kk_start = 1;
        kk_end   = p;

        gamma = opts.shifts.tmp.U(:, kk_start:kk_end)' * ...
                U_new(:, jj_start:jj_end);
        opts.shifts.tmp.RR(kk_start:kk_end, ...
                           dim_keep + jj_start:dim_keep + jj_end) = ...
            opts.shifts.tmp.RR(kk_start:kk_end, ...
                               dim_keep + jj_start:dim_keep + jj_end) + gamma;

        U_new(:, jj_start:jj_end) = U_new(:, jj_start:jj_end) - ...
            opts.shifts.tmp.U(:, kk_start:kk_end) * gamma;
    end

    % .. Gram-Schmidt to get the ortho.
    % Basis for w -> new columns in Z.
    [opts.shifts.tmp.U(:, dim_keep + jj_start:dim_keep + jj_end), ...
        opts.shifts.tmp.RR(dim_keep + jj_start:dim_keep + jj_end, ...
                           dim_keep + jj_start:dim_keep + jj_end)] = ...
        qr(U_new(:, jj_start:jj_end), 0);
end

% Update the dimensions of the subspace.
opts.shifts.tmp.dim_current        = dim_keep + dim_new;
opts.shifts.tmp.dim_total_subspace = k;
opts.shifts.tmp.dim_new            = dim_new;
end

function [eqn, opts, oper, W, Z, Y] = update_ray(eqn, opts, oper, W, Z, Y)
% Updates the Rayleigh quotients: Hproj and Eproj.
n = oper.size(eqn, opts);

if not(isfield(opts.shifts.tmp, 'Q'))
    opts.shifts.tmp.AU   = zeros(n, 0);
    opts.shifts.tmp.EU   = zeros(n, 0);
    opts.shifts.tmp.UtAU = [];
    opts.shifts.tmp.UtEU = [];
    opts.shifts.tmp.UtB  = []; % for quadratic term.
    opts.shifts.tmp.UtU  = []; % for UV' update term.
    opts.shifts.tmp.UtV  = []; % for UV' update term.
    opts.shifts.tmp.UtR  = []; % for right hand side.
    opts.shifts.tmp.Q    = [];
end

opp_type = 'N';
if eqn.type == 'N'
    opp_type = 'T';
end

if opts.shifts.naive_update_mode || isempty(opts.shifts.tmp.Q)
    % Use a naive method of updating the Rayleigh quotient by
    % recomputing everything.

    % Do the projection.
    opts.shifts.tmp.AU   = ...
        oper.mul_A(eqn, opts, opp_type, opts.shifts.tmp.U, 'N');
    opts.shifts.tmp.UtAU = opts.shifts.tmp.U' * opts.shifts.tmp.AU;

    if eqn.haveE
        opts.shifts.tmp.EU   = ...
            oper.mul_E(eqn, opts, opp_type, opts.shifts.tmp.U, 'N');
        opts.shifts.tmp.UtEU = opts.shifts.tmp.U' * opts.shifts.tmp.EU;
    end

    opts.shifts.tmp.UtB = opts.shifts.tmp.U' * eqn.BB;
    opts.shifts.tmp.UtU = opts.shifts.tmp.U' * eqn.U;
    opts.shifts.tmp.UtV = opts.shifts.tmp.U' * eqn.V;
    opts.shifts.tmp.UtR = opts.shifts.tmp.U' * W;

    if eqn.type == 'T'
        AA = opts.shifts.tmp.UtAU + ...
             opts.shifts.tmp.UtU * opts.shifts.tmp.UtV';
    else
        AA = opts.shifts.tmp.UtAU + ...
             opts.shifts.tmp.UtV * opts.shifts.tmp.UtU';
    end

    if opts.LDL_T
        opts.shifts.tmp.Hproj = ...
            full([AA, ...
                  opts.shifts.tmp.UtB * (eqn.RR \ opts.shifts.tmp.UtB'); ...
                  opts.shifts.tmp.UtR * (eqn.T * opts.shifts.tmp.UtR'), ...
                  -AA']);
    else
        opts.shifts.tmp.Hproj = ...
            full([AA, opts.shifts.tmp.UtB * opts.shifts.tmp.UtB'; ...
                  opts.shifts.tmp.UtR * opts.shifts.tmp.UtR', -AA']);
    end

    if eqn.haveE
        opts.shifts.tmp.Eproj = ...
            full(blkdiag(opts.shifts.tmp.UtEU, opts.shifts.tmp.UtEU'));
    end
else
    % Update the Rayleigh quotient by using the up/downdated
    % QR-factorization of the relevant (newest) columns of Z.
    % Currently requires matrix-vector products with A for computing
    %  A * (new columns of Q-factor) in each step.
    dim_current = opts.shifts.tmp.dim_current;
    dim_new     = opts.shifts.tmp.dim_new;

    Unew  = opts.shifts.tmp.U(:, dim_current - dim_new + 1:dim_current);

    % Compute: AUnew = A * Unew and EUnew = E * Unew;
    AUnew = oper.mul_A(eqn, opts, opp_type, Unew, 'N');

    if eqn.haveE
        EUnew = oper.mul_E(eqn, opts, opp_type, Unew, 'N');
    end

    % Update the old A*U and U'*A*U
    if size(opts.shifts.tmp.Q, 2) == 0
        opts.shifts.tmp.AU   = zeros(n, 0);
        opts.shifts.tmp.UtAU = [];
        opts.shifts.tmp.EU   = zeros(n, 0);
        opts.shifts.tmp.UtEU = [];

        if eqn.type == 'T'
            opts.shifts.tmp.UtU = Unew' * eqn.U;
        else
            opts.shifts.tmp.UtV = Unew' * eqn.V;
        end
        opts.shifts.tmp.UtB  = Unew' * eqn.BB;
    else
        opts.shifts.tmp.AU   = opts.shifts.tmp.AU * opts.shifts.tmp.Q;
        opts.shifts.tmp.UtAU = opts.shifts.tmp.Q' * ...
            (opts.shifts.tmp.UtAU * opts.shifts.tmp.Q);

        if eqn.haveE
            opts.shifts.tmp.EU   = opts.shifts.tmp.EU * ...
                opts.shifts.tmp.Q;
            opts.shifts.tmp.UtEU = opts.shifts.tmp.Q' * ...
                (opts.shifts.tmp.UtEU * opts.shifts.tmp.Q);
        end

        if eqn.type == 'T'
            opts.shifts.tmp.UtU = ...
                [opts.shifts.tmp.Q' * opts.shifts.tmp.UtU; ...
                 Unew' * eqn.U];
        else
            opts.shifts.tmp.UtV = ...
                [opts.shifts.tmp.Q' * opts.shifts.tmp.UtV; ...
                 Unew' * eqn.V];
        end
        opts.shifts.tmp.UtB  = ...
            [opts.shifts.tmp.Q' * opts.shifts.tmp.UtB; Unew' * eqn.BB];
    end

    opts.shifts.tmp.UtAU = [opts.shifts.tmp.UtAU, ...
                            opts.shifts.tmp.U(:, ...
                                              1:dim_current - dim_new)' * ...
                            AUnew; ...
                            Unew' * opts.shifts.tmp.AU, Unew' * AUnew];
    opts.shifts.tmp.AU   = [opts.shifts.tmp.AU, AUnew];

    if eqn.haveE
        opts.shifts.tmp.UtEU = [opts.shifts.tmp.UtEU, ...
                                opts.shifts.tmp.U(:, ...
                                                  1:dim_current - ...
                                                  dim_new)' * ...
                                EUnew; ...
                                Unew' * opts.shifts.tmp.EU, Unew' * EUnew];
        opts.shifts.tmp.EU   = [opts.shifts.tmp.EU, EUnew];
    end

    opts.shifts.tmp.UtR = opts.shifts.tmp.U' * W;
    if eqn.type == 'T'
        opts.shifts.tmp.UtV = opts.shifts.tmp.U' * eqn.V;
    else
        opts.shifts.tmp.UtU = opts.shifts.tmp.U' * eqn.U;
    end

    % Finally, assemble the Hamiltonian Rayleigh quotient matrix.
    if eqn.type == 'T'
        AA = opts.shifts.tmp.UtAU + ...
             opts.shifts.tmp.UtU * opts.shifts.tmp.UtV';
    else
        AA = opts.shifts.tmp.UtAU + ...
             opts.shifts.tmp.UtV * opts.shifts.tmp.UtU';
    end

    if opts.LDL_T
        opts.shifts.tmp.Hproj = ...
            [AA, ...
             opts.shifts.tmp.UtB * (eqn.RR \ opts.shifts.tmp.UtB'); ...
             opts.shifts.tmp.UtR * (eqn.T * opts.shifts.tmp.UtR'), ...
             -AA'];
    else
        opts.shifts.tmp.Hproj = ...
            [AA opts.shifts.tmp.UtB * opts.shifts.tmp.UtB'; ...
             opts.shifts.tmp.UtR * opts.shifts.tmp.UtR', -AA'];
    end

    % If needed, assemble the  Rayleigh quotient of diag(E,E') matrix.
    if eqn.haveE
        opts.shifts.tmp.Eproj = ...
            blkdiag(opts.shifts.tmp.UtEU, opts.shifts.tmp.UtEU');
    end
end
end
