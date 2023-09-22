function [Er, Ar, Br, Cr, Dr, outinfo] = mess_balanced_truncation(varargin)
% Lyapunov Balanced truncation for descriptor systems with invertible E.
%
%  [Er, Ar, Br, Cr, Dr, outinfo] = ...
%     mess_balanced_truncation(E, A, B , C, D, opts)
%
%  or
%
%  [Er, Ar, Br, Cr, Dr, outinfo] = mess_balanced_truncation(sys, opts)
%  with sys = sparss(A, B, C, D, E)
%
%  or
%
%  [Er, Ar, Br, Cr, Dr, outinfo] = mess_balanced_truncation(eqn, opts, oper)
%
% INPUTS:
%  E, A, B, C, D The mass, system, input and output matrices describing the
%                original system
%
%  sys           sparss(A, B, C, D, E)
%
%  eqn           usual eqn structure (see help mess)
%
%  oper          usual oper structure (see help mess, only non-DAE usfs allowed)
%
%  opts          options structure with substructure bt containing members
%
%      max_ord    maximum reduced order allowed
%                 (optional, defaults to size(A,1))
%
%      trunc_tol  error tolerance used for the Hankel singular value truncation
%                 (optional, defaults to 1e-5)
%
%      info       verbosity control parameter (optional):
%                 0  quiet (default)
%                 1  show iteration numbers and residuals
%                 >1 plot residual history
%                 >2 compute and show the sigma and error plots
%
%      opts can further be used to pass setting to the
%             LRADI, ADI shift computation, or the square root method
%             (optional)
%             (see corresponding routines for additional information)
%
%             opts.srm for the square root method inherits max_ord and trunc_tol
%             unless they are already given in opts.srm.
%
%             opts also has fields to control the plotting in case info>2 above:
%             opts.sigma.fmin      minimum value in the logspace for the
%                                  sigma and error plots
%             opts.sigma.fmax      maximum value in the logspace for the
%                                  sigma and error plots
%             opts.sigma.nsample   number of elements in the logspace.
%
% OUTPUTS:
% Er, Ar, Br, Cr, Dr     the reduced order model matrices
%
% outinfo.TL outinfo.TR  the left and right transformation matrices
%
% outinfo.errbound       H-infinity error upper bound
%
% outinfo.hsv            vector with the computed Hankel singular values
%

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%% check arguments
switch nargin
    case 2
        if isa(varargin{2}, 'struct')
            opts = varargin{2};
        else
            opts = struct();
            mess_err(opts, 'error_arguments', ...
                     'second argument must be a struct in case of 2 inputs');
        end
        if isa(varargin{1}, 'sparss')
            [eqn, opts, oper] = mess_wrap_sparss(varargin{1}, opts);

        else
            mess_err(opts, 'error_arguments', ...
                     'first argument must be a sparss in case of 2 inputs');
        end

    case 3
        if isa(varargin{1}, 'struct') && ...
           isa(varargin{2}, 'struct') && ...
           isa(varargin{3}, 'struct')
            eqn = varargin{1};

            opts = varargin{2};
            oper = varargin{3};
        else
            opts = struct();
            mess_err(opts, 'error_arguments', ...
                     'Each argument must be a struct in case of 3 inputs');
        end

    case 6
        if isa(varargin{6}, 'struct')
            opts = varargin{6};
        else
            opts = struct();
            mess_err(opts, 'error_arguments', ...
                     'Last argument must be a struct in case of 6 inputs');
        end

        % Problem data
        if not(issparse(varargin{1})) || not(issparse(varargin{2}))
            mess_err(opts, 'data', 'Both E and A need to be sparse.');
        end
        % needs to be saved in extra variables for octave 4.2 compatibility
        rank_1 = sprank(varargin{1});
        dimension_1 = size(varargin{1}, 1);
        if rank_1 < dimension_1
            mess_err(opts, 'data', ...
                     'Only systems with invertible E are supported at the moment');
        end

        eqn.E_ = varargin{1};
        eqn.A_ = varargin{2};
        % save non truncated matrices
        eqn.B = varargin{3};
        eqn.C = varargin{4};
        eqn.D = varargin{5};

        % operations are done by the default set of user supplied functions
        [oper, opts] = operatormanager(opts, 'default');

        % Let us avoid E if it is actually the identity.
        if norm(eqn.E_ - speye(size(eqn.E_, 1)), 'inf') == 0
            eqn.haveE = 0;
            eqn = rmfield(eqn, 'E_');
        else
            eqn.haveE = 1;
        end

    otherwise
        opts = struct();
        mess_err(opts, 'error_arguments', 'invalid number of inputs');
end

if not(isfield(opts, 'bt'))
    opts.bt = struct('tol', 1e-5, ...
                     'max_ord', oper.size(eqn, opts), ...
                     'info', 0);
end

if not(isfield(opts.bt, 'info'))
    opts.bt.info = 0;
end

%% BT tolerance and maximum order for the ROM
if not(isfield(opts, 'srm'))
    opts.srm = opts.bt;
end

if not(isfield(opts.srm, 'tol'))
    if not(isfield(opts.bt, 'tol')) || isempty(opt.bt.tol)
        opts.srm.tol = 1e-5;
    else
        opts.srm.tol = opts.bt.tol;
    end
end

if not(isfield(opts.srm, 'max_ord'))
    if not(isfield(opts.bt, 'max_ord')) || isempty(opts.bt.max_ord)
        opts.srm.max_ord = size(eqn.A_, 1);
    else
        opts.srm.max_ord = opts.bt.max_ord;
    end
end

if not(isfield(opts.srm, 'info'))
    opts.srm.info = opts.bt.info;
end

%% some control settings for the LRADI
if not(isfield(opts, 'adi'))
    % ADI tolerance and maximum iteration number
    opts.adi.maxiter = 100;
    opts.adi.rel_diff_tol = 1e-16;
    opts.adi.info = opts.bt.info;
    set_tol = true;
else
    if not(isfield(opts.adi, 'maxiter'))
        opts.adi.maxiter = 100;
    end

    if not(isfield(opts.adi, 'rel_diff_tol'))
        opts.adi.rel_diff_tol = 1e-16;
    end

    if not(isfield(opts.adi, 'info'))
        opts.adi.info = opts.bt.info;
    end

    set_tol = not(isfield(opts.adi, 'res_tol'));
end

if not(isfield(opts, 'norm'))
    opts.norm = 'fro';
end

% If not set outside, we use projection shifts
if not(isfield(opts, 'shifts')) || not(isfield(opts.shifts, 'method'))
    opts.shifts.num_desired = max(5, min(size(eqn.B, 2), size(eqn.C, 1)));
    opts.shifts.b0 = ones(oper.size(eqn, opts), 1);
    opts.shifts.method = 'projection';
end

%% Compress the RHS factors for the Lyapunov equations for robustness
% save  original matrices
B = eqn.B;
C = eqn.C;
% Truncate with machine precision error in the RHS representations, but possibly
% smaller rank when B, C are (almost) rank deficient.
eqn.B = mess_column_compression(full(eqn.B), 'N');
eqn.C = mess_column_compression(full(eqn.C), 'T');

%% Truncated controllability Gramian
eqn.type = 'N';

if set_tol % opts.adi.res_tol was not set outside

    % if users set trunc_tol == 0 we need to avoid res_tol = 0 here,
    % since otherwise LR_ADI will turn of residual checks.
    opts.adi.res_tol = min(eps / norm(eqn.B' * eqn.B), ...
                           max(eps, opts.srm.tol / 100));
end

outB = mess_lradi(eqn, opts, oper);

if outB.niter == opts.adi.maxiter
    mess_warn(opts, 'BT', ...
              ['ADI did not converge for controllability Gramian ' ...
               'factor. Reduction results may be inaccurate']);
end

if opts.bt.info > 0
    [Bm, Bn] = size(outB.Z);
    mess_fprintf(opts, 'size outB.Z: %d x %d\n', Bm, Bn);

    if opts.bt.info > 1
        plot_iter_vs_resnorm(outB.res, eqn.type, eqn.haveE);
    end
end

%% Truncated observability Gramian
eqn.type = 'T';

if set_tol % opts.adi.res_tol was not set outside

    % if users set trunc_tol == 0 we need to avoid res_tol = 0 here,
    % since otherwise LR_ADI will turn of residual checks.
    opts.adi.res_tol = min(eps / norm(eqn.C * eqn.C'), ...
                           max(eps, opts.srm.tol / 100));
end

outC = mess_lradi(eqn, opts, oper);

if outC.niter == opts.adi.maxiter
    mess_warn(opts, 'BT', ...
              ['ADI did not converge for controllability Gramian ' ...
               'factor. Reduction results may be inaccurate']);
end

if opts.bt.info > 0
    [Cm, Cn] = size(outC.Z);
    mess_fprintf(opts, 'size outC.Z: %d x %d\n', Cm, Cn);

    if opts.bt.info > 1
        plot_iter_vs_resnorm(outC.res, eqn.type, eqn.haveE);
    end
end

%% Square root method
[TL, TR, hsv] = mess_square_root_method(eqn, opts, oper, outB.Z, outC.Z);

%% compute ROM matrices
% Note that we use the original B and C since the ones in eqn have been
% truncated above.
Ar = TL' * oper.mul_A(eqn, opts, 'N', TR, 'N');
Br = TL' * B;
Cr = C * TR;
Er = eye(size(Ar, 1));

if isfield(eqn, 'D')
    Dr = eqn.D;
else
    Dr = [];
end

%% if desired, plot the approximation results
if opts.bt.info > 2
    if not(isfield(opts, 'tf_plot'))
        opts.tf_plot = struct('fmin', 1e-6, ...
                              'fmax', 1e6, ...
                              'nsample', 100, ...
                              'type', 'sigma');
    else
        if not(isfield(opts.tf_plot, 'fmin'))
            opts.tf_plot.fmin = 1e-6;
        end
        if not(isfield(opts.tf_plot, 'fmax'))
            opts.tf_plot.fmax = 1e6;
        end
        if not(isfield(opts.tf_plot, 'nsample'))
            opts.tf_plot.nsample = 100;
        end
        if not(isfield(opts.tf_plot, 'type'))
            opts.tf_plot.type = 'sigma';
        end
    end

    ROM = struct('A', Ar, 'E', Er, 'B', Br, 'C', Cr, 'D', Dr);

    % for the evaluations we need the original B and C in the eqn structure.
    eqn.B = B;
    eqn.C = C;
    mess_tf_plot(eqn, opts, oper, ROM);
end

%% construct output information
if nargout > 4
    r  = size(Ar, 1);
    nr = size(eqn.A_, 1) - length(hsv);

    outinfo = struct('TL', TL, ...
                     'TR', TR, ...
                     'errbound', 2.0 * sum(hsv(r + 1:end)) + nr * hsv(end), ...
                     'hsv', hsv);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% local function for plotting iterations vs residual norm
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_iter_vs_resnorm(out_res, eqn_type, eqn_haveE)

figure();
semilogy(out_res, 'LineWidth', 3);

if eqn_type == 'N'

    if eqn_haveE
        title('A X E^T + E X A^T = -BB^T');
    else
        title('A X + X A^T = -BB^T');
    end

elseif eqn_type == 'T'

    if eqn_haveE
        title('A^T X E + E^T X A = -C^T C');
    else
        title('A^T X + X A = -C^T C');
    end
end

xlabel('number of iterations');
ylabel('normalized residual norm');
drawnow;
end
