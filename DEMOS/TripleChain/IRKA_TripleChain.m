function IRKA_TripleChain(n1, usfs, istest)
% IRKA_TripleChain computes first order IRKA based ROM for the triple chain
%
% Usage:      IRKA_TripleChain(n1, oper, istest)
%
% Inputs:
%
% n1          length of a single chain in the model
%             (optional; defaults to 1000)
%
% usfs        the set of user supplied functions to use, or the indication to
%             call IRKA with matrices.
%             possible values: 'so_1', 'so_2', 'default', 'matrices'
%             (optional; defaults to 'so_1')
%
% istest      flag to determine whether this demo runs as a CI test or
%             interactive demo
%             (optional, defaults to 0, i.e. interactive demo)
%

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%%
narginchk(0, 3);
if nargin < 1
    n1 = 1000;
end
if nargin < 2
    usfs = 'so_1';
end
if nargin < 3
    istest = 0;
end
%%
% set usfs,  unless we want to use the matrices
opts = struct();
if strcmp(usfs, 'matrices')
    [oper, opts] = operatormanager(opts, 'so_2');
else
    [oper, opts] = operatormanager(opts, usfs);
end

% Problem data
alpha = .002;
Beta  = .002;
v     = 5;

matrices = false;
switch usfs
    case 'so_1'
        [eqn.M_, eqn.E_, eqn.K_] = ...
            triplechain_MSD(n1, alpha, Beta, v);
        s  = size(eqn.K_, 1);
        eqn.B = [zeros(s, 1); ones(s, 1)];
        eqn.C = [ones(1, s) zeros(1, s)];
    case 'so_2'
        [eqn.M_, eqn.E_, eqn.K_] = ...
            triplechain_MSD(n1, alpha, Beta, v);
        s  = size(eqn.K_, 1);
        eqn.B = [ones(s, 1); zeros(s, 1)];
        eqn.C = eqn.B';
    case 'default'
        [M_, E_, K_] = triplechain_MSD(n1, alpha, Beta, v);
        s = size(K_, 1);
        eqn.A_ = [sparse(s, s), -K_; -K_, -E_];
        eqn.E_ = [-K_, sparse(s, s); sparse(s, s), M_];
        eqn.B = [zeros(s, 1); ones(s, 1)];
        eqn.C = [ones(1, s) zeros(1, s)];
        clear M_ E_ K_ s;
    case 'matrices'
        % setup matrices for the IRKA call
        [M, E, K] = triplechain_MSD(n1, alpha, Beta, v);
        s = size(K, 1);
        B = ones(s, 1);
        Cp = B';
        Cv = zeros(1, s);
        % We will need eqn for tf_plot below
        eqn = struct('M_', M, 'E_', E, 'K_', K, ...
                     'B', [B; zeros(s, 1)], ...
                     'C', [Cp Cv]);
        matrices = true;
end
eqn.haveE = 1;
%%
opts.irka = struct('r', 30, ...
                   'maxiter', 100, ...
                   'shift_tol', 1e-3, ...
                   'h2_tol', 1e-8, ...
                   'num_prev_shifts', 5, ...
                   'num_prev_roms', 5, ...
                   'flipeig', false, ...
                   'init', 'subspace');

%%
if matrices
    [Er, Ar, Br, Cr, Dr, outinfo] = ...
        mess_tangential_irka(M, E, K, B, Cp, Cv, [], opts);
else
    [Er, Ar, Br, Cr, Dr, outinfo] = ...
        mess_tangential_irka(eqn, opts, oper);
end
ROM.E = Er;
ROM.A = Ar;
ROM.B = Br;
ROM.C = Cr;
ROM.D = Dr;

opts.tf_plot.nsample = 400;  % 400 frequency samples
opts.tf_plot.fmin = -4;      % min. frequency 1e-3
opts.tf_plot.fmax = 0;       % max. frequency 1e4
if istest
    opts.tf_plot.info = 0;
else
    opts.tf_plot.info = 3;
end
opts.tf_plot.type = 'sigma';

[out, ~, opts, ~] = mess_tf_plot(eqn, opts, oper, ROM);

if istest
    if isequal(outinfo.term_flag, 'maxiter')
        mess_err(opts, 'TEST:accuracy', ...
                 'terminated with maximum number of iterations');
    end
    if any(out.relerr > 1)
        mess_err(opts, 'TEST:accuracy', ...
                 'unexpectedly inaccurate results');

    end
end
end
