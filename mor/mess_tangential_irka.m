function [Er, Ar, Br, Cr, Dr, outinfo] = mess_tangential_irka(varargin)
% The tangential IRKA method with automatic selection of initial shifts and
% tangential directions.
%
% Call
%   [Er, Ar, Br, Cr, Dr, outinfo] = mess_tangential_irka(E, A, B, C, D, opts)
%    or
%       mess_tangential_irka(sys, opts)
%    with sys = sparss(A, B, C, D, E)
%
%   [Er, Ar, Br, Cr, Dr, outinfo] =
%       mess_tangential_irka(M, E, K, B, Cp, Cv, D, opts)
%    or
%       mess_tangential_irka(sys, opts)
%    with sys = mechss(M, E, K, B, Cp, Cv, D)
%
%   [Er, Ar, Br, Cr, Dr, outinfo] = ...
%       mess_tangential_irka(eqn, opts, oper)
%
% Inputs:
%  E, A, B, C   The mass, system, input and output matrices describing the
%               original system
%
%  M, E, K, B,  The mass, system, input and output matrices describing the
%  Cp, Cv       original system
%
%  eqn          struct contains data for equations
%
%  opts         optional options structure with substructure 'irka'
%
% Input fields in struct opts.irka:
%
%  r                 The desired reduced order (optional, default: 25)
%
%  maxiter           maximum iteration number for the IRKA iteration
%                    (optional, default: 100)
%
%  shift_tol         bound for the relative change of the IRKA shifts used
%                    as stopping criterion (optional, default: 1e-6)
%
%  h2_tol            bound for the relative change of the H2-norm compared
%                    to the last stable ROM (optional, default: 100*eps)
%
%  num_prev_shifts   number of previous shift vectors stored for cycle
%                    detection (optional, default: 5)
%
%  num_prev_ROMs     number of previous stable ROMs stored for cycle
%                    detection (optional, default: 5)
%
%  info              0  : silent
%                    1  : print status info in each IRKA step (default)
%                    >1 : compute and show the sigma and error plots
%                         (only for matrix inputs)
%
%  init              shift and direction initialization choice: (optional)
%                    'subspace' chooses a random subspace and uses it to
%                               compute projected shifts and directions
%                               from the projected EVP just like in the
%                               IRKA iteration. (default)
%                    'logspace' picks logspaced shifts in [0,1] and all
%                               ones as tangential directions
%                    'random'   picks normally distributed random shifts
%                               and tangential directions (does not fix a
%                               seed, for reproducible results the caller
%                               should do this)
%                    'rom'      an asymptotically stable initial guess for
%                               the reduced model of order r is given in
%                               opts.irka.Er, opts.irka.Ar, opts.irka.Br,
%                               opts.irka.Cr.
%
%  flipeig    flag marking whether or not to correct the signs of shifts in
%             the wrong halfplane. (optional, default: 0)
%
% Outputs:
%  Er, Ar, Br, Cr, Dr  The reduced system matrices.
%
%  outinfo             structure with members
%      S, b, c             The final shifts and tangential directions
%      TR, TL                The final truncation matrices
%      term_flag           An indicator which stopping criterion stopped IRKA
%
% NOTE: Currently only standard state space systems and descriptor systems
% with E invertible are supported, when matrices are passed in.
%
% References:
% [1] S. Gugercin, A. C. Antoulas, C. Beattie, H2 model reduction for
%     large-scale linear dynamical systems, SIAM J. Matrix Anal.
%     Appl. 30 (2) (2008) 609–638. https://doi.org/10.1137/060666123.
%

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%% decide if matrices or sparss/mechss was passed in
opts = struct();
if nargin == 2
    if isa(varargin{1}, 'sparss')
        opts = varargin{2};
        [eqn, opts, oper] = mess_wrap_sparss(varargin{1}, opts, 'default');
    elseif isa(varargin{1}, 'mechss')
        opts = varargin{2};
        [eqn, opts, oper] = mess_wrap_mechss(varargin{1}, opts, 'so_2');
    end
end

%% Choose usfs set if matrices were passed in
if nargin == 6 % default first order system
    [oper, opts] = operatormanager(opts, 'default');
elseif nargin == 8 % second order system
    [oper, opts] = operatormanager(opts, 'so_2');
end

%% Fill equation structure if matrices were passed in
if nargin == 6
    eqn.A_ = varargin{2};
    if isempty(varargin{1})
        eqn.E_ = speye(size(varargin{2}, 1));
    else
        eqn.E_ = varargin{1};
        eqn.haveE = true;
    end
    eqn.B = full(varargin{3});
    eqn.C = full(varargin{4});
    eqn.D = full(varargin{5});
    opts = varargin{6};
elseif nargin == 8
    eqn.M_ = varargin{1};
    eqn.E_ = varargin{2};
    eqn.K_ = varargin{3};
    eqn.haveE = true;
    eqn.B = [full(varargin{4}); zeros(size(varargin{4}))];
    eqn.C = [full(varargin{5}), full(varargin{6})];
    eqn.D = full(varargin{7});
    opts = varargin{8};
end

%% handle D term
D = [];

if nargin == 3
    eqn = varargin{1};
    opts = varargin{2};
    oper = varargin{3};
    if isfield(eqn, 'D')
        D = eqn.D;
    end
end
Dr = D;

%% check field opts.irka
if not(isfield(opts, 'irka')) || not(isstruct(opts.irka))
    opts.irka = [];
end

if notfield_or_empty(opts.irka, 'r')
    opts.irka.r = 25;
end

if notfield_or_empty(opts.irka, 'num_prev_shifts')
    opts.irka.num_prev_shifts = 5;
end

if notfield_or_empty(opts.irka, 'num_prev_ROMs')
    opts.irka.num_prev_ROMs = 5;
end

if notfield_or_empty(opts.irka, 'maxiter')
    opts.irka.maxiter = 100;
end

if notfield_or_empty(opts.irka, 'init')
    opts.irka.init = 'subspace';
end

if notfield_or_empty(opts.irka, 'shift_tol')
    opts.irka.shift_tol = 1e-6;
end

if notfield_or_empty(opts.irka, 'h2_tol')
    opts.irka.h2_tol = 100 * eps;
end

if notfield_or_empty(opts.irka, 'info')
    opts.irka.info = 1;
end

if notfield_or_empty(opts.irka, 'flipeig')
    opts.irka.flipeig = false;
end

%% Initialize used usfs
[result, eqn, opts, oper] = oper.init(eqn, opts, oper, 'A', 'E');

if not(result)
    mess_err(opts, 'data', 'Equation  data seems to be incomplete');
end

[eqn, opts, oper] = oper.mul_A_pre(eqn, opts, oper);
[eqn, opts, oper] = oper.mul_E_pre(eqn, opts, oper);
[eqn, opts, oper] = oper.sol_ApE_pre(eqn, opts, oper);

%% Initialization
n = oper.size(eqn, opts);
m = size(eqn.B, 2);
p = size(eqn.C, 1);

% preallocation of space for storing previous iterates and shifts
S_old = zeros(opts.irka.r, opts.irka.num_prev_shifts);
Er_old = cell(opts.irka.num_prev_ROMs, 1);
Ar_old = cell(opts.irka.num_prev_ROMs, 1);
Br_old = cell(opts.irka.num_prev_ROMs, 1);
Cr_old = cell(opts.irka.num_prev_ROMs, 1);

% new field saving whether the reduced order model is stable
isstab = zeros(opts.irka.maxiter, 1);

% initialize termination flag
outinfo.term_flag = [];

% Choose the initial shifts and tangential directions
initial_rom = false;

switch opts.irka.init

    case 'subspace'
        U = get_initial_subspace(n, opts.irka.r);
        [X, S, Y] = eig(full(U' * (oper.mul_A(eqn, opts, 'N', U, 'N'))), ...
                        full(U' * (oper.mul_E(eqn, opts, 'N', U, 'N'))));

        M = diag(ones(opts.irka.r, 1) ./ sqrt(diag(Y' * X)));
        Y = Y * M;
        X = X * M;
        b = (Y' * (U' * eqn.B)).';
        c = (eqn.C * U) * X;

        [S, perm] = mess_make_proper(diag(S));
        b = b(:, perm);
        c = c(:, perm);

    case 'random'
        S = abs(randn(opts.irka.r, 1));
        b = randn(m, opts.irka.r);
        c = randn(p, opts.irka.r);

    case 'logspace'
        S = logspace(0, 1, opts.irka.r);
        b = ones(m, opts.irka.r);
        c = ones(p, opts.irka.r);

    case 'rom'
        Er = opts.irka.Er;
        Ar = opts.irka.Ar;
        Br = opts.irka.Br;
        Cr = opts.irka.Cr;
        [X, S, Y] = eig(Ar, Er);
        if any(real(diag(S)) >= 0)
            mess_err(opts, 'unstable', ...
                     ['The initial guess for the reduced system must ', ...
                      'be asymptotically stable!']);
        end

        M = diag(ones(opts.irka.r, 1) ./ sqrt(diag(Y' * X)));
        Y = Y * M;
        X = X * M;
        b = (Y.' * Br).';
        c = Cr * X;

        [S, perm] = mess_make_proper(diag(S));
        b = b(:, perm);
        c = c(:, perm);

        initial_rom = true;
end

%% Start iteration
for iter = 1:opts.irka.maxiter

    %% save previous shifts
    if iter <= opts.irka.num_prev_shifts
        S_old(:, iter) = S;
    else
        S_old(:, mod(iter, opts.irka.num_prev_shifts) + 1) = S;
    end

    %% save previous ROM into the buffer of recently seen stable ROMs
    if (initial_rom && (iter == 1)) || ((iter > 1) && isstab(iter - 1))

        if initial_rom && (iter == 1)
            old = 1;
        end

        if (iter > 1) && isstab(iter - 1)
            num_stab_ROMs = sum(isstab);
            if num_stab_ROMs <= opts.irka.num_prev_ROMs
                old = num_stab_ROMs;
            else
                old = mod(num_stab_ROMs, opts.irka.num_prev_ROMs) + 1;
            end
        end

        Er_old{old} = Er;
        Ar_old{old} = Ar;
        Br_old{old} = Br;
        Cr_old{old} = Cr;
    end

    %% Compute projection subspaces
    TR = zeros(n, opts.irka.r);
    TL = zeros(n, opts.irka.r);
    k = 1;

    while k < (opts.irka.r + 1)
        x = oper.sol_ApE(eqn, opts, 'N', S(k), 'N', eqn.B * b(:, k), 'N');
        y = oper.sol_ApE(eqn, opts, 'T', S(k), 'T', eqn.C' * c(:, k), 'N');

        if not(imag(S(k)) == 0)
            TR(:, k:k + 1) = [real(x), imag(x)];
            TL(:, k:k + 1) = [real(y), imag(y)];
            k = k + 2;
        else
            TR(:, k) = real(x);
            TL(:, k) = real(y);
            k = k + 1;
        end
    end

    % find orthonormal bases
    [TR, ~] = qr(TR, 0);
    [TL, ~] = qr(TL, 0);

    %% Biorthogonalize TR,TL in the E inner product
    Er = (TL' * oper.mul_E(eqn, opts, 'N', TR, 'N'));
    [U, Sigma, Q] = svd(Er);
    Sigma = diag(ones(opts.irka.r, 1) ./ sqrt(diag(Sigma)));
    TL = TL * U * Sigma;
    TR = TR * Q * Sigma;

    %% Compute new ROM
    Ar = TL' * oper.mul_A(eqn, opts, 'N', TR, 'N');
    Br = TL' * eqn.B;
    Cr = eqn.C * TR;
    Er = eye(opts.irka.r); % by construction

    %% Update interpolation points/tangential directions
    % compute eigendecomposition
    [X, S, Y] = eig(Ar);
    S = diag(S);

    % ensure the correct scaling of left and right eigenvectors
    M = diag(ones(opts.irka.r, 1) ./ sqrt(diag(Y' * X)));
    Y = Y * M;
    X = X * M;

    % make sure all shifts are in the correct half plane.
    wrongsign = find(real(S) > 0);

    if not(isempty(wrongsign))

        if opts.irka.flipeig
            mess_warn(opts, 'unstable', ...
                      ['IRKA step %d : %d non-stable reduced eigenvalues ' ...
                       'have been flipped.\n'], iter, length(wrongsign));
        else
            mess_warn(opts, 'unstable', ...
                      ['IRKA step %d : %d non-stable reduced eigenvalues ' ...
                       'detected.\n'], iter, length(wrongsign));
        end
    else
        isstab(iter) = 1;
    end

    if opts.irka.flipeig
        S(wrongsign) = -S(wrongsign);
    end

    % update tangential directions
    b = (Y.' * Br).';
    c = Cr * X;

    % make sure complex conjugate shifts come in consecutive pairs such
    % that real basis extension works properly
    [S, perm] = mess_make_proper(S);
    b = b(:, perm);
    c = c(:, perm);

    %% compute convergence indices
    % maximum pointwise relative change of the shifts for the last
    % num_prev_shifts shift vectors
    shiftchg = realmax;

    for shift_iter = 1:min(opts.irka.num_prev_shifts, iter)

        shiftchg_iter = norm((S - S_old(:, shift_iter)) ./ S, 'inf');

        if shiftchg_iter < shiftchg
            shiftchg = shiftchg_iter;
        end
    end

    % relative H2-change
    % Computation of the change only makes sense if we have an old model
    % and both, the current and old, models are stable
    if (initial_rom && isstab(iter)) || ...
       ((iter > 1) && (isstab(iter) && any(isstab(1:iter - 1))))
        romchg = realmax;

        for rom_iter = 1:min(opts.irka.num_prev_ROMs, sum(isstab) - 1)

            romchg_iter = mess_h2_rom_change(Er, Ar, Br, Cr, ...
                                             Er_old{rom_iter}, ...
                                             Ar_old{rom_iter}, ...
                                             Br_old{rom_iter}, ...
                                             Cr_old{rom_iter}, 1);
            if romchg_iter < romchg
                romchg = romchg_iter;
            end
        end
    else
        if iter == 1
            romchg = 1.0;
        else
            romchg = Inf;
        end
    end

    %% If desired print status message
    if opts.irka.info
        mess_fprintf(opts, ...
                     ['IRKA step %3d, rel. chg. shifts = %e , rel. H2-norm', ...
                      ' chg. ROM = %e\n'], ...
                     iter, shiftchg, romchg);
    end

    %% evaluate stopping criteria
    if shiftchg < opts.irka.shift_tol

        outinfo.term_flag = 'shift_tol';

        if opts.irka.info
            mess_fprintf(opts, ['IRKA terminated due to relative change ', ...
                                'of shifts criterion.\n\n']);
        end

        break
    end

    if romchg < opts.irka.h2_tol

        outinfo.term_flag = 'h2_tol';

        if opts.irka.info
            mess_fprintf(opts, ['IRKA terminated due to relative change ', ...
                                'of ROMs criterion.\n\n']);
        end

        break
    end
end

if (iter == opts.irka.maxiter) && isempty(outinfo.term_flag)
    mess_warn(opts, 'convergence', ...
              'IRKA: No convergence in %d iterations.\n', opts.irka.maxiter);
    outinfo.term_flag = 'maxiter';
end

if opts.irka.info > 1

    ROM = struct('A', Ar, 'E', Er, 'B', Br, 'C', Cr, 'D', D);
    if not(isfield(opts, 'tf_plot'))
        opts.tf_plot = struct();
    end
    if not(isfield(opts.tf_plot, 'fmin'))
        opts.tf_plot.fmin = -6;
    end
    if not(isfield(opts.tf_plot, 'fmax'))
        opts.tf_plot.fmax = 6;
    end
    if not(isfield(opts.tf_plot, 'nsample'))
        opts.tf_plot.nsample = 100;
    end
    if not(isfield(opts.tf_plot, 'info'))
        opts.tf_plot.info = opts.irka.info;
    end
    if not(isfield(opts.tf_plot, 'type'))
        opts.tf_plot.type = 'sigma';
    end

    [~, eqn, opts, oper] = mess_tf_plot(eqn, opts, oper, ROM);
end

[eqn, opts, oper] = oper.mul_A_post(eqn, opts, oper);
[eqn, opts, oper] = oper.mul_E_post(eqn, opts, oper);

oper.sol_ApE_post(eqn, opts, oper);

outinfo.S = S;
outinfo.b = b;
outinfo.c = c;

outinfo.TR = TR;
outinfo.TL = TL;

end

function bool = notfield_or_empty(mystruct, myfield)
bool = not(isfield(mystruct, myfield)) || isempty(mystruct.(myfield));
end
