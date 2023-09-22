function [eqn, opts, oper] = sol_ApE_pre_so_iter(eqn, opts, oper)
% function [eqn, opts, oper] = sol_ApE_pre_so_iter(eqn, opts, oper)
% To simplify matters in sol_ApE we add a field eqn.E_ holding the
% identity matrix when we do not have an E matrix already.

% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)

if not(eqn.haveE)
    eqn = fake_E_default(eqn);
end

%% Check input Parameters
if not(isfield(eqn, 'haveE'))
    eqn.haveE = false;
end
if eqn.haveE
    if not(isfield(eqn, 'E_')) || not(isfield(eqn, 'K_')) || not(isfield(eqn, 'M_'))
        mess_err(opts, 'error_arguments', 'Field eqn.M_ or eqn.E_ or eqn.K_ is not defined');
    end
else
    if not(isfield(eqn, 'M_'))
        mess_err(opts, 'error_arguments', 'Field eqn.M_ is not defined');
    end
    if not(isfield(eqn, 'K_'))
        mess_err(opts, 'error_arguments', 'Field eqn.K_ is not defined');
    end
    if not(isfield(eqn, 'E_'))
        mess_err(opts, 'error_arguments', 'Field eqn.E_ is not defined');
    end
end

% Check for the solver used
if isfield(opts.usfs.so_iter, 'method_ApE')
    if not(exist(opts.usfs.so_iter.method_ApE, 'file') == 2)
        mess_err(opts, 'control_data', ['iterative solver method field ''method_ApE''', ...
                                        ' is an unsupported solver.']);
    end
else
    mess_warn(opts, 'control_data', ['iterative solver method field ''method_ApE''', ...
                                     ' is unset. Falling back to GMRES.']);
    opts.usfs.so_iter.method_ApE = 'gmres';
end

%  Restart size for GMRES
if strcmpi(opts.usfs.so_iter.method_ApE, 'gmres')

    if isfield(opts.usfs.so_iter, 'restIter')
        if opts.usfs.so_iter.restIter < 0
            mess_err(opts, 'control_data', ['GMRES restart iterations value', ...
                                            'is invalid']);
        end
    else
        mess_warn(opts, 'control_data', ['GMRES restart iterations not', ...
                                         ' found. Falling back to default']);
        opts.usfs.so_iter.restIter = 25;
    end
end

%% Pre-defined preconditioner

if not(isfield(opts.usfs.so_iter, 'PApE_R'))

    if not(isfield(opts.usfs.so_iter, 'PApE_L'))

        mess_warn(opts, 'control_data', ['No preconditioner for ApE could be found.', ...
                                         ' Switching to ILU']);
        form_ApE = @(alpha, p)ApEfun(alpha, p, eqn);
        [L, U] = ilu(form_ApE(opts.usfs.so_iter.alpha, opts.usfs.so_iter.p_));
        opts.usfs.so_iter.PApE_L = L;
        opts.usfs.so_iter.PApE_R = U;
    else

        opts.usfs.so_iter.PApE_R = [];

    end

end

end

function Y = ApEfun(alpha, p, eqn)
% Y = ApEfun(alpha, p, eqn);
% This is a function handler that accepts the parameter alpha and p, to return
% the matrix ApE_ for solving a second order system re-shaping it as a first
% order system.
Y = [-alpha * eqn.K_ + p * eqn.K_, eqn.K_ - alpha * eqn.E_ + p * alpha * eqn.M_; ...
     -eqn.K_ + p * alpha * eqn.M_, -eqn.E_ + alpha * eqn.E_ + p * eqn.M_];
end
