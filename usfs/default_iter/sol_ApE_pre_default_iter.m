function [eqn, opts, oper] = sol_ApE_pre_default_iter(eqn, opts, oper)
% function [eqn, opts, oper] = sol_ApE_pre_default_iter(eqn, opts, oper)
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
    if not(isfield(eqn, 'E_')) || not(isfield(eqn, 'A_'))
        mess_err(opts, 'error_arguments', 'Field eqn.E_ or eqn.A_ is not defined');
    end
else
    if not(isfield(eqn, 'A_'))
        mess_err(opts, 'error_arguments', 'Field eqn.A_ is not defined');
    end
    if isfield(eqn, 'E_')
        mess_err(opts, 'equation_data', ['Detected eqn.E_ where eqn.haveE ' ...
                                         'is 0. You need to set haveE = true or delete E_.']);
    end
end

% Check for the solver used
if isfield(opts.usfs.default_iter, 'method_ApE')
    if not(exist(opts.usfs.default_iter.method_ApE, 'file') == 2)
        mess_err(opts, 'control_data', ['iterative solver method field ''method_ApE''', ...
                                        ' is an unsupported solver.']);
    end
else
    mess_warn(opts, 'control_data', ['iterative solver method field ''method_ApE''', ...
                                     ' is unset. Falling back to GMRES.']);
    opts.usfs.default_iter.method_ApE = 'gmres';
end

%  Restart size for GMRES
if strcmpi(opts.usfs.default_iter.method_ApE, 'gmres')

    if isfield(opts.usfs.default_iter, 'restIter')
        if opts.usfs.default_iter.restIter < 0
            mess_err(opts, 'control_data', ['GMRES restart iterations value', ...
                                            'is invalid']);
        end
    else
        mess_warn(opts, 'control_data', ['GMRES restart iterations not', ...
                                         ' found. Falling back to default']);
        opts.usfs.default_iter.restIter = 25;
    end
end

%% Pre-defined preconditioner

if not(isfield(opts.usfs.default_iter, 'PApE_R'))

    if not(isfield(opts.usfs.default_iter, 'PApE_L'))

        mess_warn(opts, 'control_data', ['No preconditioner for ApE could be found.', ...
                                         ' Switching to ILU']);
        [L, U] = ilu(-eqn.A_ + opts.usfs.default_iter.p_ * eqn.E_);
        opts.usfs.default_iter.PApE_L = L;
        opts.usfs.default_iter.PApE_R = U;
    else

        opts.usfs.default_iter.PApE_R = [];

    end

end

end
