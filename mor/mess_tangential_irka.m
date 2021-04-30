function [Er, Ar, Br, Cr, S, b, c, V, W, term_flag] = mess_tangential_irka(varargin)
% The tangential IRKA method with automatic selection of initial shifts and
% tangential directions.
%
% Call
%   [Er,Ar,Br,Cr,S,b,c,V,W,term_flag] = mess_tangential_irka(E,A,B,C,opts)
%                                       or
%                                     = mess_tangential_irka(sys,opts)
%                                       with sys = sparss(A,B,C,[],E)
%
%   [Er,Ar,Br,Cr,S,c,V,W,term_flag] = mess_tangential_irka(M,E,K,B,Cp,Cv,opts)
%                                       or
%                                   = mess_tangential_irka(sys,opts)
%                                     with sys = mechss(M,E,K,B,Cp,Cv,[])
%
%   [Er,Ar,Br,Cr,S,c,V,W,term_flag] = mess_tangential_irka(eqn,opts,oper)
%
% Inputs:
%  E,A,B,C    The mass, system, input and output matrices describing the
%             original system
%
%  M,E,K,B,   The mass, system, input and output matrices describing the
%  Cp,Cv      original system
%
%   eqn       struct contains data for equations
%
%  opts       optional options structure with substructure 'irka'
%
% Input fields in struct opts.irka:
%
%  r                 The desired reduced order (optional, default: 25)
%
%  maxiter           maximum iteration number for the IRKA iteration
%                    (optional, default: 25)
%
%  shift_tol         bound for the relative change of the IRKA shifts used as
%                    stopping criterion (optional, default: 1e-2)
%
%  h2_tol            bound for the relative change of the H2-norm compared to
%                    the last stable ROM (optional, default: 100*eps)
%
%  num_prev_shifts   number of previous shift vectors stored for cycle
%                    detection (optional, default: 5)
%
%  num_prev_ROMs     number of previous stable ROMs stored for cycle
%                    detection (optional, default: 5)
%
%  info              0  : silent (default)
%                    1  : print status info in each IRKA step
%                    >1 : compute and show the sigma and error plots
%                         (only for matrix inputs)
%
%  init              shift and direction initialization choice: (optional)
%                    'subspace' chooses a random subspace and uses it to compute
%                               projected shifts and directions from the
%                               projected EVP just like in the irka iteration.
%                               (default)
%                    'logspace' picks logspaced shifts in [0,1] and all ones
%                               as tangential directions
%                    'random'   picks normally distributed random shifts and
%                               tangential directions
%                    'rom'      an asymptotically stable initial guess for the
%                               reduced model of order r is given in
%                               opts.irka.Er, opts.irka.Ar, opts.irka.Br,
%                               opts.irka.Cr.
%
%  flipeig    flag marking whether or not to correct the signs of shifts in
%             the wrong halfplane. (optional, default: 1)
%
% Outputs:
%  Er,Ar,Br,Cr   The reduced system matrices.
%  S,b,c         The final shifts and tangential directions
%  V,W           The final truncation matrices
%  term_flag     An indicator which stopping criterion stopped IRKA
%
% NOTE: Currently only standard state space systems and descriptor systems
% with E invertible are supported.
%
% References:
% [1] S. Gugercin, A. C. Antoulas, C. Beattie, H2 model reduction for
%     large-scale linear dynamical systems, SIAM J. Matrix Anal. Appl. 30 (2)
%     (2008) 609–638. https://doi.org/10.1137/060666123.
%

%
% This file is part of the M-M.E.S.S. project 
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2021 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%% decide if matrices or sparss/mechss was passed in
if nargin == 2
    if isa(varargin{1}, 'sparss')
        [eqn, oper] = mess_wrap_sparss(varargin{1}, 'default');
        opts = varargin{2};
    elseif isa(varargin{1}, 'mechss')
        [eqn, oper] = mess_wrap_mechss(varargin{1}, 'so_2');
        opts = varargin{2};
    end
end

%% Choose usfs set if matrices were passed in
if nargin == 5 % default first order system
    oper = operatormanager('default');
elseif nargin == 7 % second order system
    oper = operatormanager('so_2');
end

%% Fill equation structure if matrices were passed in
if nargin == 5
    eqn.A_ = varargin{2};
    if isempty(varargin{1})
        eqn.E_ = speye(size(varargin{2}, 1));
    else
        eqn.E_ = varargin{1};
        eqn.haveE = 1;
    end
    eqn.B = full(varargin{3});
    eqn.C = full(varargin{4});
    opts = varargin{5};
elseif nargin == 7
    eqn.M_ = varargin{1};
    eqn.E_ = varargin{2};
    eqn.K_ = varargin{3};
    eqn.haveE = 1;
    eqn.B = [full(varargin{4}); zeros(size(varargin{4}))];
    eqn.C = [full(varargin{5}), full(varargin{6})];
    opts = varargin{7};
end
D = [];

%%
if nargin == 3
    eqn = varargin{1};
    opts = varargin{2};
    oper = varargin{3};
    if isfield(eqn, 'D')
        D = eqn.D;
    end
end

%% check field opts.irka
if not(isfield(opts, 'irka')) || not(isstruct(opts.irka))
    opts.irka = [];
end

if not(isfield(opts.irka, 'r')) || isempty(opts.irka.r)
    opts.irka.r = 25;
end

if not(isfield(opts.irka, 'num_prev_shifts')) || ...
        isempty(opts.irka.num_prev_shifts)
    opts.irka.num_prev_shifts = 5;
end
if not(isfield(opts.irka, 'num_prev_ROMs')) || ...
        isempty(opts.irka.num_prev_ROMs)
    opts.irka.num_prev_ROMs = 5;
end
if not(isfield(opts.irka, 'maxiter')) || isempty(opts.irka.maxiter)
    opts.irka.maxiter = 100;
end
if not(isfield(opts.irka, 'init')) || isempty(opts.irka.init)
    opts.irka.init = 'subspace';
end
if not(isfield(opts.irka, 'shift_tol')) || isempty(opts.irka.shift_tol)
    opts.irka.shift_tol = 1e-6;
end
if not(isfield(opts.irka, 'h2_tol')) || isempty(opts.irka.h2_tol)
    opts.irka.h2_tol = 100 * eps;
end
if not(isfield(opts.irka, 'info')) || isempty(opts.irka.info)
    opts.irka.info = 1;
end
if not(isfield(opts.irka, 'flipeig')) || isempty(opts.irka.flipeig)
    opts.irka.flipeig = 0;
end

%% Initialize used usfs
[result, eqn, opts, oper] = oper.init(eqn, opts, oper, 'A', 'E');
if not(result)
    error('MESS:data', 'Equation  data seems to be incomplete');
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
term_flag = [];

% Choose the initial shifts and tangential dirctions
initial_rom = 0;
switch opts.irka.init
    case 'subspace'
        U = get_initial_subspace(n, opts.irka.r);
        [T, S] = eig(full(U' * (oper.mul_A(eqn, opts, 'N', U, 'N'))), ...
            full(U' * (oper.mul_E(eqn, opts, 'N', U, 'N'))));
        [S, perm] = mess_make_proper(diag(S));
        S = S';
        T = T(:, perm);
        b = (T \ U' * eqn.B).';
        c = eqn.C * U * T;
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
        [T, S] = eig(Ar, Er);
        if any(real(diag(S)) >= 0)
            error(['The initial guess for the reduced system must ', ...
                'be asymptotically stable!']);
        end
        [S, perm] = mess_make_proper(diag(S));
        S = conj(S);
        T = T(:, perm);
        b = (T \ Br).';
        c = Cr * T;
        initial_rom = 1;
end

%% Start iteration
for iter = 1:opts.irka.maxiter

    %% save previous shifts
    if iter <= opts.irka.num_prev_shifts
        S_old(:, iter) = S;
    else
        S_old(:, mod(iter, opts.irka.num_prev_shifts)+1) = S;
    end

    %% save previous ROM
    if (initial_rom && iter == 1) || (iter > 1 && isstab(iter - 1))
        if initial_rom && iter == 1
            old = 1;
        end
        if iter > 1 && isstab(iter-1)
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
    V = zeros(n, opts.irka.r);
    W = zeros(n, opts.irka.r);
    j = 1;
    while (j < opts.irka.r+1)
        x = oper.sol_ApE(eqn, opts, 'N', S(j), 'N', eqn.B*b(:, j), 'N');
        y = oper.sol_ApE(eqn, opts, 'T', S(j), 'T', eqn.C'*c(:, j), 'N');
        if (imag(S(j)) ~= 0)
            V(:, j:j+1) = [real(x), imag(x)];
            W(:, j:j+1) = [real(y), imag(y)];
            j = j + 2;
        else
            V(:, j) = real(x);
            W(:, j) = real(y);
            j = j + 1;
        end
    end
    % find orthonormal bases
    [V, ~] = qr(V, 0);
    [W, ~] = qr(W, 0);

    %% Biorthogonalize V,W in the E inner product
    Er = (W' * oper.mul_E(eqn, opts, 'N', V, 'N'));
    [U, Sigma, Q] = svd(Er);
    Sigma = diag(ones(opts.irka.r, 1)./sqrt(diag(Sigma)));
    W = W * U * Sigma;
    V = V * Q * Sigma;

    %% Compute new ROM

    Ar = W' * oper.mul_A(eqn, opts, 'N', V, 'N');
    Br = W' * eqn.B;
    Cr = eqn.C * V;
    Er = eye(opts.irka.r); % by construction

    %% Update interpolation points/tangential directions
    [T, S] = eig(Ar);
    [S, perm] = mess_make_proper(diag(S));
    T = T(:, perm);
    % make sure all shifts are in the coorect half plane.
    wrongsign = find(real(S) > 0);
    if not(isempty(wrongsign))
        if (opts.irka.flipeig)
            warning('MESS:IRKA:unstable', ...
                ['IRKA step %d : %d non-stable reduced eigenvalues have ', ...
                'been flipped.\n'], iter, length(wrongsign));
        else
            warning('MESS:IRKA:unstable', ...
                'IRKA step %d : %d non-stable reduced eigenvalues detected.\n', ...
                iter, length(wrongsign));
        end
    else
        isstab(iter) = 1;
    end
    if (opts.irka.flipeig)
        S(wrongsign) = -S(wrongsign);
    end
    % update tangential directions
    b = (T \ Br).';
    c = Cr * T;

    %% compute convergence indices
    % maximum pointwise relative change of the shifts for the last num_prev_shifts
    % shift vectors
    shiftchg = realmax;
    for shift_iter = 1:min(opts.irka.num_prev_shifts, iter)
        shiftchg_iter = norm((S-S_old(:, shift_iter))./S, 'inf');
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
        for rom_iter = 1:min(opts.irka.num_prev_ROMs, sum(isstab)-1)
            romchg_iter = mess_h2_rom_change(Er, Ar, Br, Cr, ...
                Er_old{rom_iter}, Ar_old{rom_iter}, Br_old{rom_iter}, ...
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
        fprintf(['IRKA step %3d, rel. chg. shifts = %e , rel. H2-norm ', ...
            'chg. ROM = %e\n'], iter, shiftchg, romchg);
    end

    %% evaluate stopping criteria
    if shiftchg < opts.irka.shift_tol
        term_flag = 'shift_tol';
        if opts.irka.info
            fprintf(['IRKA terminated due to relative change ', ...
                'of shifts criterion.\n\n']);
        end
        break;
    end
    if romchg < opts.irka.h2_tol
        term_flag = 'h2_tol';
        if opts.irka.info
            fprintf(['IRKA terminated due to relative change ', ...
                'of ROMs criterion.\n\n']);
        end
        break;
    end
end
if ((iter == opts.irka.maxiter) && isempty(term_flag))
    warning('MESS:IRKA:convergence', ...
        'IRKA: No convergence in %d iterations.\n', opts.irka.maxiter);
    term_flag = 'maxiter';
end

if opts.irka.info > 1 && nargin > 3
    ROM = struct('A', Ar, 'E', Er, 'B', Br, 'C', Cr, 'D', D);
    if not(isfield(opts, 'sigma')), opts.sigma = struct(); end
    if not(isfield(opts.sigma, 'fmin')), opts.sigma.fmin = -6; end
    if not(isfield(opts.sigma, 'fmax')), opts.sigma.fmax = 6; end
    if not(isfield(opts.sigma, 'nsample')), opts.sigma.nsample = 100; end
    if not(isfield(opts.sigma, 'info')), opts.sigma.info = opts.irka.info; end
    [~, eqn, opts, oper] = mess_sigma_plot(eqn, opts, oper, ROM);
end

[eqn, opts, oper] = oper.mul_A_post(eqn, opts, oper);
[eqn, opts, oper] = oper.mul_E_post(eqn, opts, oper);
oper.sol_ApE_post(eqn, opts, oper);
