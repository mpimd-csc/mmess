function KSM_Rail(k, type, transformed, eqtype, SISO, space, istest)
% KSM_RAIL solves Lyapunov or Riccati equations for the revised
% Oberwolfach steel profile benchmark problem.
%
% INPUTS:
% k            refinement level for the rail data. See
%              mess_get_linear_rail for allowed values.
%              (optional, default: 1)
%
% type         selects Lyapunov ('LE') or Riccati ('CARE') equation
%
% transformed  number of degrees of freedom per spatial
%              direction, i.e. problem size is n^2
%              (optional, default: 50)
%
% eqntype      selects primal ('N') or dual ('T') equation
%              (optional, default: 'N')
%
% space        the type of Krylov subspace to use.
%                  'EK'  extended Krylov
%                  'RK'  rational Krylov
%              (optional, default: 'RK')
%
% istest       only used on our continuous integration
%              infrastructure. Adds another accuracy check when true.
%              (optional, default: 'false')
%
% OUTPUTS:     none
%

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

if nargin < 7
    istest = false;
end

if nargin < 6
    space = 'EK';
end

if nargin < 1
    k = 1;
end
% fetch rail model of appropriate size
eqn = mess_get_linear_rail(k);
eqn.haveE = true;

if nargin < 5
    SISO = false;
end
% For the SISO test pick one column of B and one row of C
if SISO
    eqn.B = sum(eqn.B, 2);
    eqn.C = sum(eqn.C, 1);
end

%%
if nargin < 3
    transformed = false;
end
if transformed
    % Note that this is for demonstration purposes only and should never be
    % done explicitly, otherwise.
    eqn.A_ = eqn.E_ \ eqn.A_;
    eqn.B = eqn.E_ \ eqn.B;
    eqn.haveE = false;
    eqn.E_ = speye(size(eqn.A_, 1));
end

%%
% set operation
% the Krylov projection method require E = I so we currently make
% them require the state_space_transformed_default usfs, which
% work on a transformed version of the system without forming the
% transformed matrices.
opts = struct();
[oper, opts] = operatormanager(opts, 'state_space_transformed_default');

%%
if nargin < 4
    eqtype = 'N';
end
eqn.type = eqtype;
opts.KSM.maxiter = 100;
opts.KSM.res_tol = 1e-6;
opts.KSM.comp_res = 1;
opts.KSM.space = space;
switch space
    case 'RK'
        opts.KSM.type_shifts = 'real';
end

if nargin < 2
    type = 'CARE';
elseif strcmp(type, 'LE') && strcmp(type, 'CARE')
    mess_err(opts, 'illegal_input', ...
             'type must be either ''LE'' or ''CARE''');
end
opts.KSM.type_eqn = type;
opts.KSM.CARE_shifts = 'Ritz';
opts.KSM.symmetric = true;
opts.KSM.trunc_tol = 1e-17;
opts.KSM.trunc_info = 1;
opts.KSM.explicit_proj = false;

opts.KSM.info = 1;

opts.norm = 'fro';

mess_fprintf(opts, '\n');
mess_fprintf(opts, 'Computing solution factors\n');

t_KSM = tic;
[out, eqn, opts, ~] = mess_KSM(eqn, opts, oper);
toc(t_KSM);
%%
% this shouldn't be computed for large problems
if k > 4
    return
end

X = out.Z * out.D * out.Z';
mess_fprintf(opts, '\n');  % \n at the beginning of a string breaks
% spellchecking

mess_fprintf(opts, ...
             ['The final internal and real normalized residual norms ', ...
              'after %d iterations are:\n'], out.niter);

if strcmp(eqn.type, 'T')
    if  strcmp(opts.KSM.type_eqn, 'CARE')
        mess_fprintf(opts, '\n');
        nrm = norm(eqn.A_' * X + ...
                   X * eqn.A_ - ...
                   (X * eqn.B) * (eqn.B' * X) + ...
                   eqn.C' * eqn.C, 'fro') / ...
              norm(eqn.C' * eqn.C, 'fro');
        mess_fprintf(opts, ['CARE Frobenius norm: %10g (real) \t %10g ', ...
                            '(internal)\n'], ...
                     nrm, out.res(end));

        gnrm = norm(eqn.A_' * X * eqn.E_ + ...
                    eqn.E_' * X * eqn.A_ - ...
                    eqn.E_' * (X * eqn.B) * (eqn.B' * X) * eqn.E_ + ...
                    eqn.C' * eqn.C, 'fro') / ...
               norm(eqn.C' * eqn.C, 'fro');
        mess_fprintf(opts, ['gCARE Frobenius norm: %10g (real) \t %10g ', ...
                            '(internal)\n'], ...
                     gnrm, out.res(end));
    else
        mess_fprintf(opts, '\n');
        nrm = norm(eqn.A_' * X + X * eqn.A_ + eqn.C' * eqn.C, 'fro') / ...
              norm(eqn.C' * eqn.C, 'fro');
        mess_fprintf(opts, ['LE Frobenius norm: %10g (real) \t %10g ', ...
                            '(internal)\n'], ...
                     nrm, out.res(end));

        gnrm = norm(eqn.A_' * X * eqn.E_ + ...
                    eqn.E_' * X * eqn.A_ + ...
                    eqn.C' * eqn.C, 'fro') / ...
               norm(eqn.C' * eqn.C, 'fro');
        mess_fprintf(opts, ['gLE Frobenius norm: %10g (real) \t %10g ', ...
                            '(internal)\n'], ...
                     gnrm, out.res(end));
    end
else
    if  strcmp(opts.KSM.type_eqn, 'CARE')
        mess_fprintf(opts, '\n');
        nrm = norm(eqn.A_ * X + ...
                   X * eqn.A_' - ...
                   (X * eqn.C') * (eqn.C * X) + ...
                   eqn.B * eqn.B', 'fro') / ...
              norm(eqn.B * eqn.B', 'fro');
        mess_fprintf(opts, ['CARE Frobenius norm: %10g (real) \t %10g ', ...
                            '(internal)\n'], ...
                     nrm, out.res(end));

        gnrm = norm(eqn.A_ * X * eqn.E_' + ...
                    eqn.E_ * X * eqn.A_' - ...
                    eqn.E_ * (X * eqn.C') * (eqn.C * X) * eqn.E_' + ...
                    eqn.B * eqn.B', 'fro') / ...
               norm(eqn.B * eqn.B', 'fro');
        mess_fprintf(opts, ['gCARE Frobenius norm: %10g (real) \t %10g ', ...
                            '(internal)\n'], ...
                     gnrm, out.res(end));
    else
        mess_fprintf(opts, '\n');
        nrm = norm(eqn.A_ * X + X * eqn.A_' + eqn.B * eqn.B', 'fro') / ...
              norm(eqn.B * eqn.B', 'fro');
        mess_fprintf(opts, ['LE Frobenius norm: %10g (real) \t %10g ' ...
                            '(internal)\n'], ...
                     nrm, out.res(end));

        gnrm = norm(eqn.A_ * X * eqn.E_' + ...
                    eqn.E_ * X * eqn.A_' + ...
                    eqn.B * eqn.B', 'fro') / ...
               norm(eqn.B * eqn.B', 'fro');
        mess_fprintf(opts, ...
                     'gLE Frobenius norm: %10g (real) \t %10g (internal)\n', ...
                     gnrm, out.res(end));
    end
end

if istest
    if not(transformed)
        if (abs(gnrm - out.res(end)) / gnrm) > 1e-4
            mess_err(opts, 'failure', 'test failed');
        end
    else
        if (abs(nrm - out.res(end)) / nrm) > 1e-4
            mess_err(opts, 'failure', 'test failed');
        end
    end
end
