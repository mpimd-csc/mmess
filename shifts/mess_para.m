function [p, out, eqn, opts, oper] = mess_para(eqn, opts, oper)
%
%  Estimation of suboptimal ADI shift parameters for the matrix (operator) F=A
%
%  Calling sequence:
%
%    [p, out, eqn, opts, oper] = mess_para(eqn, opts, oper)
%
%  Input:
%
%    eqn       structure contains data A, E, B, C, K
%
%    opts      struct contains parameters for the algorithm
%
%    oper      contains function handles with operations for A and E
%
%  Output:
%
%    p         an opts.shifts.num_desired- or opts.shifts.num_desired+1-vector of
%              suboptimal ADI parameters;
%
%    out       outputstructure potentially containing the following fields
%              (depending on the method used): 
%    out.err_code  Error code = 1, if Ritz values with positive real parts
%                  have been encountered; otherwise, err_code = 0;
%    out.rw        vector containing the Ritz values;
%    out.Hp        Hessenberg matrix in Arnoldi process w.r.t. F;
%    out.Hm        Hessenberg matrix in Arnoldi process w.r.t. inv(F);
%    out.Vp        Orthogonal matrix in Arnoldi process w.r.t. F;
%    out.Vm        Orthogonal matrix in Arnoldi process w.r.t. inv(F);
%
%    eqn, opts, oper  potentially altered versions of the above inputs.
%
% Input fields in struct eqn:
%
%   eqn.B       dense (n x m1) matrix B
%
%   eqn.C       dense (m2 x n) matrix C
%
%   eqn.U       dense (n x m3) matrix U
%               (required if eqn.V is present)
%
%   eqn.V       dense (n x m3) matrix V
%               (required if eqn.U is present)
%
%   eqn.type    possible  values: 'N', 'T'
%               determining whether (N) or (T) is solved
%               (optional)
%
%   eqn.haveE   possible  values: 0, 1, false, true
%               if haveE = 0: matrix E is assumed to be the identity
%               (optional)
%
%   Depending on the operator chosen by the operatormanager, additional
%   fields may be needed. For the "default", e.g., eqn.A_ and eqn.E_ hold
%   the A and E matrices. For the second order types these are given
%   implicitly by the M, D, K matrices stored in eqn.M_, eqn.E_ and eqn.K_,
%   respectively.
%
% Input fields in struct opts:
%   opts.shifts.num_desired   possible  values: integer > 0
%                             number of shifts that should be computed
%                             2*num_desired < num_Ritz + num_hRitz is required
%                             (optional, default: 25)
%
%   opts.shifts.num_Ritz      possible  values: integer > 0
%                             number of Arnoldi steps w.r.t. F for
%                             heuristic shift computation
%                             num_Ritz < n is required
%                             (optional, default: 50)
%
%   opts.shifts.num_hRitz     possible  values: integer > 0
%                             number of Arnoldi steps w.r.t. inv(F) for
%                             heuristic shift computation
%                             num_hRitz < n is required
%                             (optional, default: 25)
%
%   opts.shifts.b0            (n x 1) array
%                             start vector for Arnoldi algorithm for
%                             heuristic shift computation
%                             (optional, default: ones(n, 1))
%
%   opts.shifts.info          possible  values: 0, 1, false, true
%                             turn output of used shifts before the first
%                             iteration step on (1) or off (0)
%                             (optional, default: 0)
%
%   opts.shifts.method        possible  values:
%                             'heuristic', ('heur', 'Penzl', 'penzl')
%                                for Penzl's heuristics.
%                             'wachspress', ('Wachspress')
%                                for asymptotically optimal Wachspress
%                                selection.
%                             'projection'
%                                for adaptively updated projection shifts.
%                             method for shift computation
%                             in case of 'projection' new shifts are
%                             computed during the iteration steps,
%                             otherwise the shifts are reused cyclically
%                             (optional, default: 'heuristic')
%
%   opts.shifts.truncate      possible values: scalar >= 1.0
%                             truncation tolerance to drop exceptionally large
%                             and small Ritz values (e.g. used in second order
%                             cases, where (p^2*M + p*E + K) may otherwise
%                             numerically loose the information about either
%                             M or K in finite precision). Ritz values larger
%                             than the given value and smaller than its
%                             reciprocal are truncated.
%                             (optional, default: [])
%
%   opts.shifts.banned        array of shift parameter values that the
%                             ADI will not use.
%                             shift parameter computation will remove all
%                             shifts in a neighborhood of the banned
%                             shifts.
%                             use opts.shifts.banned_tol as relative
%                             tolerance for the neighborhood size
%                             (optional, default: [])
%
%   opts.shifts.banned_tol    possible  values: scalar >= 0
%                             relative tolerance for the neighborhood
%                             size around banned shifts
%                             (optional, default: 1e-4)
%
%   opts.shifts.implicitVtAV  possible values: true, false
%                             decides whether A*V is reconstructed
%                             implicitly in 'projection' method, unused
%                             otherwise.
%                             (optional, default: true)
%
%  Remarks:
%
%    Typical values are opts.shifts.num_desired = 10..40,
%    opts.shifts.num_Ritz = 20..80, opts.shifts.num_hRitz = 10..40.
%    The harder the problem is the large values are necessary.
%    Larger values mostly result in a faster convergence, but also in a
%    larger memory requirement.
%    However, for "well-conditioned" problems small values of
%    opts.shifts.num_desired can lead to the optimal performance.
%    In case of the projection shifts, a natural selection for l0
%    is the number of columns, i.e. normally the rank, of the right hand side.
%
%  References:
%
%  [1] T. Penzl.
%      LYAPACK (Users' Guide - Version 1.0).
%      1999.
%
%  [2] P. Kürschner, Efficient low-rank solution of large-scale matrix
%      equations, Dissertation, Otto-von-Guericke-Universität, Magdeburg,
%      Germany, shaker Verlag, ISBN 978-3-8440-4385-3 (Apr. 2016).
%      URL http://hdl.handle.net/11858/00-001M-0000-0029-CE18-2
%
%   uses operator function size
%   and indirectly requires
%   (heuristic and wachspress shifts)
%    - size, sol_A, mul_A, sol_E, mul_E in mess_arn
%   (projection shifts)
%    - mul_A, mul_E  in mess_projection_shifts (projection shifts

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

% Input data not completely checked!

%% check data

if not(isfield(opts.shifts,'method'))
    opts.shifts.method='heuristic';
    warning('MESS:control_data', ...
        ['Missing shift parameter selection method. ', ...
         'Switching to default: heuristic shifts']);
end

if not(isfield(opts,'shifts')) || not(isstruct(opts.shifts))
    warning('MESS:control_data',...
        ['shift parameter control structure missing. ', ...
         'Switching to defaults: ', ...
         'num_desired = 25, num_Ritz = 50, num_hRitz = 25.']);
    opts.shifts.num_desired = 25;
    opts.shifts.num_Ritz = 50;
    opts.shifts.num_hRitz = 25;
else
    if not(isfield(opts.shifts,'num_desired')) || ...
            not(isnumeric(opts.shifts.num_desired))
        warning('MESS:control_data',...
            ['Missing or Corrupted opts.shifts.num_desired field.', ...
            'Switching to default: 25']);
        opts.shifts.num_desired = 25;
    end
    if strcmp(opts.shifts.method,'heur') && ...
       (not(isfield(opts.shifts,'num_Ritz')) || ...
       not(isnumeric(opts.shifts.num_Ritz)))
        warning('MESS:control_data',...
            ['Missing or Corrupted opts.shifts.num_Ritz field.', ...
            'Switching to default: 50']);
        opts.shifts.num_Ritz = 50;
    end
    if strcmp(opts.shifts.method,'heur') && ...
       (not(isfield(opts.shifts,'num_hRitz')) || ...
        not(isnumeric(opts.shifts.num_hRitz)))
        warning('MESS:control_data',...
            ['Missing or Corrupted opts.shifts.num_hRitz field.', ...
            'Switching to default: 25']);
        opts.shifts.num_hRitz = 25;
    end
end

if not(isfield(eqn, 'haveE')), eqn.haveE = 0; end

if not(isfield(eqn, 'type')), eqn.type = 'N'; end

[result, eqn, opts, oper] = oper.init(eqn, opts, oper, 'A','E');

if not(result)
    error('MESS:control_data', ...
        'system data is not completely defined or corrupted');
end

out.err_code = 0;

rosenbrock = 0;
if isfield(opts,'rosenbrock') && isstruct(opts.rosenbrock) && ...
        isfield(opts.rosenbrock,'tau')
    rosenbrock = 1;
    if opts.rosenbrock.stage == 1
        pc = -1 / (2 * opts.rosenbrock.tau);
        taugamma = 1;
    else % p = 2
        taugamma = (opts.rosenbrock.tau * opts.rosenbrock.gamma);
        pc = ( - 0.5) / taugamma;
    end
end

bdf = 0;
if isfield(opts,'bdf') && isstruct(opts.bdf) && ...
        isfield(opts.bdf, 'tau') && isfield(opts.bdf, 'beta')
    bdf = 1;
    pc = -1 / (2 * opts.bdf.tau * opts.bdf.beta);
end


if not(isfield(opts.shifts, 'banned')) ...
        || not(isnumeric(opts.shifts.banned))
    opts.shifts.banned = [];
elseif not(isfield(opts.shifts, 'banned_tol')) ...
        || not(isnumeric(opts.shifts.banned_tol)) ...
        || not(isscalar(opts.shifts.banned_tol))
    opts.shifts.banned_tol = 1e-4;
end

if not(isfield(opts.shifts,'recursion_level')) ...
        || not(isnumeric(opts.shifts.recursion_level)) ...
        || not(isscalar(opts.shifts.recursion_level))
    opts.shifts.recursion_level = 0;
end

%% initialize usfs
[eqn,opts, oper] = oper.mul_A_pre(eqn, opts, oper);
[eqn,opts, oper] = oper.mul_E_pre(eqn, opts, oper);
[eqn,opts, oper] = oper.sol_A_pre(eqn, opts, oper);
[eqn,opts, oper] = oper.sol_E_pre(eqn, opts, oper);

%%
switch opts.shifts.method
    case {'heur', 'heuristic', 'penzl', 'Penzl'}
        %%
        if isfield(oper, 'get_ritz_vals')
            [rw, out.Hp, out.Hm, out.Vp, out.Vm, eqn, opts, oper] = ...
                oper.get_ritz_vals(eqn, opts, oper);
        else
            [rw, out.Hp, out.Hm, out.Vp, out.Vm, eqn, opts, oper] = ...
                mess_get_ritz_vals(eqn, opts, oper);
        end

        p = mess_mnmx(rw, opts.shifts.num_desired);

    case {'wachspress', 'Wachspress'}
        %%
        if isfield(oper, 'get_ritz_vals')
            if nargout < 4
                rw = oper.get_ritz_vals(eqn, opts, oper);
            else
                [rw, out.Hp, out.Hm, out.Vp, out.Vm, eqn, opts, oper] = ...
                    oper.get_ritz_vals(eqn, opts, oper);
            end
        else
            if nargout < 4
                rw = mess_get_ritz_vals(eqn, opts, oper);
            else
                [rw, out.Hp, out.Hm, out.Vp, out.Vm, eqn, opts, oper] = ...
                    mess_get_ritz_vals(eqn, opts, oper);
            end
        end

        a = min(abs(real(rw)));
        b = max(abs(real(rw)));
        alpha = atan(max(imag(rw) ./ real(rw)));

        if not(isfield(opts.shifts, 'wachspress'))
            opts.shifts.wachspress = 'T';
        end
        switch opts.shifts.wachspress
            case 'N'
                p = mess_wachspress_n(a,b,alpha,opts.shifts.num_desired);
            case 'T'
                if isfield(opts,'nm') && ...
                        isfield(opts.nm,'inexact') &&...
                        isa(opts.nm.inexact,'char')
                    tol = opts.adi.outer_tol;
                else
                    tol=opts.adi.res_tol;
                end
                p = mess_wachspress(a, b, alpha, tol);
            otherwise
                error('MESS:shift_method',...
                      'wachspress selector needs to be either ''T'' or ''N''');
        end
    case 'projection'
        if isfield(eqn, 'G')
            U = eqn.G;
        elseif eqn.type == 'N'
            U = eqn.B;
        else
            U = eqn.C';
        end
        if issparse(U)
            U = full(U);
        end
        p = [];
        i = 1;
        while isempty(p)
            if bdf
                AU = (opts.bdf.tau * opts.bdf.beta) * ...
                    oper.mul_ApE(eqn, opts, eqn.type, pc, eqn.type, U, 'N');

            elseif rosenbrock
                AU = taugamma * ...
                    oper.mul_ApE(eqn, opts, eqn.type, pc, eqn.type, U, 'N');
                    
                if isfield(eqn, 'haveUV') && eqn.haveUV
                    if eqn.type == 'N'
                        AU = AU + eqn.U * (eqn.V' * U);
                    else
                        AU = AU + eqn.V * (eqn.U' * U);
                    end
                end
                
            else
                AU = oper.mul_A(eqn, opts, eqn.type, U, 'N');

                if isfield(eqn, 'haveUV') && eqn.haveUV
                    if eqn.type == 'T'
                        AU = AU + eqn.V * (eqn.U' * U);
                    else
                        AU = AU + eqn.U * (eqn.V' * U);
                    end
                end
            end

            if isfield(oper,'get_ritz_vals')
                p = oper.get_ritz_vals(eqn, opts, oper, U, AU, []);
            else
                p = mess_projection_shifts(eqn, opts, oper, U, AU, []);
            end

            if isempty(p)
                if  (i < 5)
                    warning('MESS:mess_para', ...
                        ['Could not compute initial projection shifts. ',...
                        'Going to retry with random right hand side.']);
                    U = rand(size(U));
                else
                    error('MESS:mess_para', ...
                        'Could not compute initial projection shifts.');
                end
            end
            
            i = i + 1;
        end
    otherwise
        error('MESS:shift_method', ...
            'unknown shift computation method requested.');
end

%% check computed shifts
% check for banned shifts
for j = 1 : length(opts.shifts.banned)
    critical_shifts = abs(p - opts.shifts.banned(j)) ...
        < opts.shifts.banned_tol * max(abs(p));
    p(critical_shifts) = p(critical_shifts) ...
        - opts.shifts.banned_tol * 2;
    %     p = p(not(critical_shifts));
end
if isempty(p) % if all shifts banned try again with double amount
    if opts.shifts.recursion_level < 2
        warning('MESS:mess_para', 'All computed shifts are banned. Retrying');
        num_desired = opts.shifts.num_desired;
        opts.shifts.num_desired = num_desired * 2;
        opts.shifts.recursion_level = opts.shifts.recursion_level + 1;
        [p , ~, eqn, opts, oper] = mess_para(eqn, opts, oper);
        opts.shifts.num_desired = num_desired;
        opts.shifts.recursion_level = opts.shifts.recursion_level - 1;
    end
end

p = mess_make_proper(p);

%% finalize usfs
[eqn, opts, oper] = oper.mul_A_post(eqn, opts, oper);
[eqn, opts, oper] = oper.mul_E_post(eqn, opts, oper);
[eqn, opts, oper] = oper.sol_A_post(eqn, opts, oper);
[eqn, opts, oper] = oper.sol_E_post(eqn, opts, oper);
end
