function KSM_FDM(type, n, eqntype, space, istest)
% KSM_FDM solves Lyapunov or Riccati equations for a FDM discretized
% convective heat transport problem.
%
% INPUTS:
% type      selects Lyapunov ('LE') or Riccati ('CARE') equation
%
% n         number of degrees of freedom per spatial
%           direction, i.e. problem size is n^2
%           (optional, default: 50)
%
% eqntype   selects primal ('N') or dual ('T') equation
%           (optional, default: 'N')
%
% space     the type of Krylov subspace to use.
%              'EK'  extended Krylov
%              'RK'  rational Krylov
%           (optional, default: 'RK')
%
% istest    only used on our continuous integration
%           infrastructure. Adds another accuracy check when true.
%           (optional, default: 'false')
%
% OUTPUTS:  none
%

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

if nargin < 5
    istest = false;
end

if nargin < 4
    space = 'RK';
end

if nargin < 3
    type = 'CARE';
end

if nargin < 2
    n = 50;
end

if nargin < 1
    eqntype = 'N';
end

% Problem data
% 2d convective heat equation
eqn.A_ = fdm_2d_matrix(n, '10*x', '100*y', '0');
eqn.B = fdm_2d_vector(n, '.1<x<=.3');
eqn.C = fdm_2d_vector(n, '.7<x<=.9');
eqn.C = eqn.C';
eqn.haveE = false;
%%
% set operation
% the Krylov projection method require E = I so we currently make
% them require the state_space_transformed_default usfs, which
% work on a transformed version of the system without forming the
% transformed matrices. (No transformation is required for this model)
opts = struct();
[oper, opts] = operatormanager(opts, 'state_space_transformed_default');
%%

eqn.type = eqntype;

opts.LDL_T = false;
opts.norm = 'fro';
opts.shifts.info = 1;
opts.KSM.maxiter = 100;
opts.KSM.space = space;

% rational Krylov needs some settings for the shifts
switch space
    case 'RK'
        opts.KSM.type_shifts = 'complex';
        opts.KSM.init_shifts(1) = norm(eqn.A_, 1);
        opts.KSM.init_shifts(2) = opts.KSM.init_shifts(1) / condest(eqn.A_);
        opts.KSM.CARE_shifts = 'Ritz';
end

% set solver for the projected problems depending on the
% availability of the control toolbox (MATLAB) / package (OCTAVE).
switch upper(type)
    case 'CARE'
        if exist('icare', 'file') % this is TMW's recommended solver
            opts.KSM.projection.meth = 'icare';
        elseif exist('care', 'file') % use care if MATLAB is too
            % old for icare ore we are on OCTAVE (with control package)
            opts.KSM.projection.meth = 'care';
        else  % if neither is available (i.e. we are on plain
            % MATLAB or OCTAVE) use dense Kleinman-Newton shipped
            % with M-M.E.S.S. as a fall back.
            opts.KSM.projection.meth = 'mess_dense_nm';
        end

    case 'LE'
        if exist('lyap', 'file') % If we have a control toolbox /
            % control package use lyap
            opts.KSM.projection.meth = 'lyap';
        else  % otherwise fall back to our shipped 2Solve method.
            opts.KSM.projection.meth = 'lyap2solve';
        end

end

opts.KSM.type_eqn = type;
opts.KSM.symmetric = false;
opts.KSM.trunc_tol = 1e-16;
opts.KSM.explicit_proj = false;
opts.KSM.res_tol = 1e-6;
opts.KSM.comp_res = 1;

opts.KSM.info = 1;

fprintf('\nComputing solution factorization.\n');

solve = tic;
[out, eqn, opts, ~] = mess_KSM(eqn, opts, oper);
toc(solve);

%%
% Let's compare the estimated residual with the real one.
% This shouldn't be computed for large problems
if n <= 80
    X = out.Z * out.D * out.Z';

    if strcmp(opts.KSM.type_eqn, 'LE')
        if strcmp(eqn.type, 'N')
            nrm = norm(eqn.A_ * X + X * eqn.A_' + ...
                       eqn.B * eqn.B', 'fro') / ...
                  norm(eqn.B' * eqn.B, 'fro');
        elseif strcmp(eqn.type, 'T')
            nrm = norm(eqn.A_' * X + X * eqn.A_ + ...
                       eqn.C' * eqn.C, 'fro') / ...
                  norm(eqn.C * eqn.C', 'fro');
        end
    elseif strcmp(opts.KSM.type_eqn, 'CARE')
        if strcmp(eqn.type, 'N')
            nrm = norm(eqn.A_ * X + X * eqn.A_' - ...
                       (X * eqn.C') * (eqn.C * X) + ...
                       eqn.B * eqn.B', 'fro') / ...
                  norm(eqn.B' * eqn.B, 'fro');
        else
            nrm = norm(eqn.A_' * X + X * eqn.A_ - ...
                       (X * eqn.B) * (eqn.B' * X) + ...
                       eqn.C' * eqn.C, 'fro') / ...
                  norm(eqn.C * eqn.C', 'fro');
        end
    end

    mess_fprintf(opts, ...
                 ['The final internal and real normalized residual ', ...
                  'norms after %3d iterations are:\n %10e (int)\t ', ...
                  '%10e (real)\n'], ...
                 out.niter, out.res(end), nrm);

end

if istest
    if out.res(end) > opts.KSM.res_tol || nrm > opts.KSM.res_tol
        mess_err(opts, 'failure', 'test failed');
    end
end
