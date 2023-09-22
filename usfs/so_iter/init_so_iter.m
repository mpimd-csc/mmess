function [result, eqn, opts, oper] = init_so_iter(eqn, opts, oper, flag1, flag2)

% function [result, eqn, opts, oper] = init_so_iter(eqn, opts, oper, flag1, flag2)
%
% The function returns true or false if data for A_ and E_
% resp. flag1 and flag2  are available and correct in structure
% eqn.
%
%   Input:
%
%   eqn             structure with data
%   opts            structure containing parameter for the algorithm
%   oper            structure contains function handles for operation with A and E
%   flag1           'A'/'E' to check if A or E is in eqn
%   flag2           'A'/'E' to check if A or E is in eqn
%
%   Output:
%
%   result          1 if data corresponding to flag1 (and flag2)
%   is available , 0 if data is not available
%   eqn             structure with data
%   opts            structure containing parameter for the algorithm
%   oper            structure contains function handles for operation with A and E
%
%   This function does not use other default functions.
%
%   This function calls two other functions checkA and checkE
%   implemented at the end.
%
%   The function checkA(eqn) proves if a field 'A_' is included in
%   the structure eqn and if the field 'A_' is numeric and
%   quadratic.
%
%   The function checkE(eqn) proves if a field 'E_' is included in
%   the structure eqn and if the field 'E_' is numeric and
%   quadratic.
%   If the structure does not include a field E, a new field 'E_'
%   is defined as a sparse identity matrix by size of field 'A_'.
%

% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%% Calculate Alpha
mfun = @(x) (eqn.M_ * x + eqn.E_ * (eqn.K_ \ (eqn.E_ * x)) / 4);
colD = size(eqn.E_, 1);
opts.usfs.so_iter.alpha = 1 / (2 * eigs(mfun, colD, eqn.E_, 1, 'LR'));

%% Check input Parameters
% Check if M is positive definite

if isfield(eqn, 'M_')
    if issymmetric(eqn.M_)
        opts.usfs.so_iter.M_is_sym = true;
        if not(any(diag(eqn.M_ <= 0)))
            opts.usfs.so_iter.M_is_spd = true;
        else
            opts.usfs.so_iter.M_is_spd = false;
        end
    else
        opts.usfs.so_iter.M_is_sym = false;
        opts.usfs.so_iter.M_is_spd = false;
    end
end

% Check if E is positive definite

if isfield(eqn, 'E_')
    if issymmetric(eqn.E_)
        opts.usfs.so_iter.E_is_sym = true;
        if not(any(diag(eqn.E_ <= 0)))
            opts.usfs.so_iter.E_is_spd = true;
        else
            opts.usfs.so_iter.E_is_spd = false;
        end
    else
        opts.usfs.so_iter.E_is_sym = false;
        opts.usfs.so_iter.E_is_spd = false;
    end
end

% Check if K is positive definite

if isfield(eqn, 'K_')
    if issymmetric(eqn.K_)
        opts.usfs.so_iter.K_is_sym = true;
        if not(any(diag(eqn.K_ <= 0)))
            opts.usfs.so_iter.K_is_spd = true;
        else
            opts.usfs.so_iter.K_is_spd = false;
        end
    else
        opts.usfs.so_iter.K_is_sym = false;
        opts.usfs.so_iter.K_is_spd = false;
    end
end

% Verification checks, to determine if a proper alpha has been selected.

% E is positive definite iff:
% i) M > 0, which is true by assumption,
% ii) K–α2M > 0,
% iii) Since E is built symmetric plus the previous 2 conditions, it becomes spd
if not(any(diag(eqn.M_ <= 0)))
    if not(any(diag((eqn.K_ - 2 * opts.usfs.so_iter.alpha * eqn.M_) <= 0)))
        opts.usfs.so_iter.E_is_spd = true;
    else
        opts.usfs.so_iter.E_is_spd = false;
    end
else
    opts.usfs.so_iter.E_is_spd = false;
end

% A is negative definite iff
% i) –αK < 0, which is true by assumption,
% ii) – D + α ( M + 1/4 (D*(K^-1)* D ) < 0 ,
% iii) A is built NO symmetric plus the previous 2 conditions, it becomes NO spd.
if any(diag((-opts.usfs.so_iter.alpha * eqn.K_) < 0))
    if any(-eqn.E_ + opts.usfs.so_iter.alpha * (eqn.M_ + ...
                                                (eqn.E_ * (eqn.K_^-1) * eqn.E_) / 4))
        opts.usfs.so_iter.A_is_spd = 0;
    else
        opts.usfs.so_iter.A_is_spd = -1;
        mess_warn(opts, 'control_data', 'A is not negative definite');
    end
else
    opts.usfs.so_iter.A_is_spd = -1;
end

% Required fields for the iterative solver
if not(isfield(opts, 'usfs')) || not(isfield(opts.usfs, 'so_iter'))

    mess_warn(opts, 'control_data', [' The ''so_iter'' usfs need the ', ...
                                     '''opts.usfs.so_iter'' substructure to be present.']);

end

% Check for the solver used
if isfield(opts.usfs.so_iter, 'method_A')
    if not(exist(opts.usfs.so_iter.method_A, 'file') == 2)
        mess_err(opts, 'control_data', ['Iterative solver method field ''method_A''', ...
                                        ' is an unsupported solver.']);
    end
else
    mess_warn(opts, 'control_data', ['Iterative solver method field ''method_A''', ...
                                     ' is unset. Falling back to GMRES.']);
    opts.usfs.so_iter.method_A = 'gmres';

end

if isfield(opts.usfs.so_iter, 'method_E')
    if not(exist(opts.usfs.so_iter.method_E, 'file') == 2)
        mess_err(opts, 'control_data', ['Iterative solver method field ''method_E''', ...
                                        ' is an unsupported solver.']);
    end
else
    if opts.usfs.so_iter.K_is_spd && opts.usfs.so_iter.M_is_spd
        mess_warn(opts, 'control_data', ['Iterative solver method field ''method_E''', ...
                                         ' is unset. Falling back to PCG.']);
        opts.usfs.so_iter.method_E = 'pcg';

    else
        mess_warn(opts, 'control_data', ['Iterative solver method field ''method_E''', ...
                                         ' is unset. Falling back to GMRES.']);
        opts.usfs.so_iter.method_E = 'gmres';
    end
end

% Required residual tolerance for stopping the iterative solver
if isfield(opts.usfs.so_iter, 'res_tol')
    if opts.usfs.so_iter.res_tol < 0
        mess_err(opts, 'control_data', ['Iterative solver residual tolerance value', ...
                                        'is invalid']);
    end
else
    mess_warn(opts, 'control_data', ['Iterative solver residual tolerance value not', ...
                                     ' found. Falling back to default']);
    opts.usfs.so_iter.res_tol = 1e-12;

end

% Number of iterations allowed

if isfield(opts.usfs.so_iter, 'max_iter')
    if opts.usfs.so_iter.max_iter < 0
        mess_err(opts, 'control_data', ['Iterative solver max. iterations', ...
                                        'is invalid']);
    end
else
    mess_warn(opts, 'control_data', ['Iterative solver max. iterations not', ...
                                     ' found. Falling back to default']);
    opts.usfs.so_iter.max_iter = oper.size(eqn, opts, oper);

end

% Restart size for GMRES
if strcmpi(opts.usfs.so_iter.method_A, 'gmres') || ...
        strcmpi(opts.usfs.so_iter.method_E, 'gmres')

    if isfield(opts.usfs.so_iter, 'restIter')
        if opts.usfs.so_iter.restIter < 0
            mess_err(opts, 'control_data', ['GMRES restart iterations value', ...
                                            'is invalid']);
        end
    else
        mess_warn(opts, 'control_data', ['GMRES restart iterations value not', ...
                                         ' found. Falling back to default']);
        opts.usfs.so_iter.restIter = 25;
    end

end

%% Operations
na = nargin;
if isfield(eqn, 'LTV')
    mess_warn(opts, 'not_implemented', ['''so_iter'' does not yet support', ...
                                        ' LTV systems.']);
end
if na < 4
    mess_err(opts, 'control_data', 'Number of input Arguments are at least 4');

elseif nargin == 4     % result = init_so_iter(eqn, flag1);
    switch flag1
        case {'A', 'a'}
            [eqn, opts, result] = checkA(eqn, opts);

        case {'E', 'e'}
            [eqn, opts, result] = checkE(eqn, opts);

        otherwise
            mess_err(opts, 'control_data', 'flag1 has to be ''A_'' or ''E_''');
    end

elseif nargin == 5     % result = init_so_iter(eqn,flag1,flag2);
    switch flag1
        case {'A', 'a'}
            [eqn, opts, result] = checkA(eqn, opts);

            switch flag2
                case {'A', 'a'}
                case {'E', 'e'}
                    [eqn, opts, resultE] = checkE(eqn, opts);
                    result = result && resultE;

                otherwise
                    mess_err(opts, 'control_data', 'flag2 has to be ''A'' or ''E''');
            end

        case {'E', 'e'}
            [eqn, opts, result] = checkE(eqn, opts);

            switch flag2
                case {'A', 'a'}
                    [eqn, opts, resultA] = checkA(eqn, opts);
                    result = result && resultA;

                case {'E', 'e'}

                otherwise
                    mess_err(opts, 'control_data', 'flag2 has to be ''A'' or ''E''');

            end

        otherwise
            mess_err(opts, 'control_data', 'flag1 has to be ''A'' or ''E''');

    end

end
end

% Check data for A_
function [eqn, opts, result] = checkA(eqn, opts)
% This function returns the changed structure eqn and a boolean
% value result (1 if 'A_' is in structure eqn and a numeric, symmetric and
% quadratic field, 0 otherwise).
% This function also defines the preconditioner for A_ by using
% ILU function.

result = isfield(eqn, 'M_') && isfield(eqn, 'E_') && isfield(eqn, 'K_');
if result
    if not(opts.usfs.so_iter.M_is_spd)
        mess_err(opts, 'control_data', ['The set of so_iter functions has only been' ...
                                        'implemented for systems where M is symmetric and positive definite']);
    end
    if not(opts.usfs.so_iter.E_is_spd)
        mess_err(opts, 'control_data', ['The set of so_iter functions has only been' ...
                                        'implemented for systems where E is symmetric and positive definite']);
    end
    if not(opts.usfs.so_iter.K_is_spd)
        mess_err(opts, 'control_data', ['The set of so_iter functions has only been' ...
                                        'implemented for systems where K is symmetric and positive definite']);
    end
end

if not(isfield(opts.usfs.so_iter, 'PA_R'))

    if not(isfield(opts.usfs.so_iter, 'PA_L'))

        mess_warn(opts, 'control_data', ['No preconditioner for A could be found.', ...
                                         ' Switching to ICHOL/ILU']);
        form_A = @(alpha)afun(alpha, eqn);
        [L, U] = ilu(form_A(opts.usfs.so_iter.alpha));
        opts.usfs.so_iter.PA_L = L;
        opts.usfs.so_iter.PA_R = U;

    else

        opts.usfs.so_iter.PA_R = [];

    end

end

end

% Check data for E_
function [eqn, opts, result] = checkE(eqn, opts)
% This function returns the changed structure eqn and a boolean
% value result (1 if 'E_' is in structure eqn and a numeric, symmetric and
% quadratic field, 0 otherwise).
% This function also defines the preconditioner for E_ by using ICHOL or
% ILU function, depending on whether E_ is symmetric and positive
% definite, or not.

if not(isfield(eqn, 'haveE'))
    eqn.haveE = false;
end % Does this still apply?

% if not(eqn.haveE)
%     if isfield(eqn, 'E_')
%         mess_err(opts,'equation_data', ['Detected eqn.E_ where eqn.haveE ' ...
%                             'is 0. You need to set haveE = true or ' ...
%                             'delete E_.']);
%     else
%         result = true;
%     end
% else

result = isfield(eqn, 'M_') && isfield(eqn, 'K_');
if result
    if not(opts.usfs.so_iter.M_is_spd)
        mess_err(opts, 'control_data', ['The set of so_iter functions has only been' ...
                                        'implemented for systems where M is symmetric and positive definite']);
    end
    if not(opts.usfs.so_iter.K_is_spd)
        mess_err(opts, 'control_data', ['The set of so_iter functions has only been' ...
                                        'implemented for systems where K is symmetric and positive definite']);
    end
end

% Preconditioner for E
if not(isfield(opts.usfs.so_iter, 'PE_R'))

    if not(isfield(opts.usfs.so_iter, 'PE_L'))

        mess_warn(opts, 'control_data', ['No preconditioner for E could be', ...
                                         ' found. Switching to ICHOL/ILU']);
        form_E = @(alpha)efun(alpha, eqn);
        if result
            S = ichol(form_E(opts.usfs.so_iter.alpha));
            opts.usfs.so_iter.PE_L = S;
            opts.usfs.so_iter.PE_R = S';

        else
            [L, U] = ilu(form_E(opts.usfs.so_iter.alpha));
            opts.usfs.so_iter.PE_L = L;
            opts.usfs.so_iter.PE_R = U;
        end
    else

        opts.usfs.so_iter.PE_R = [];

    end

end

end

function Y = afun(alpha, eqn)
% Y = afun(alpha, eqn)
% This is a function handler that accepts the parameter alpha and returns
% the matrix A_ for solving a second order system re-shaping it as a first
% order system.
Y = [-alpha * eqn.K_, eqn.K_ - alpha * eqn.E_; -eqn.K_, -eqn.E_ + alpha * eqn.E_];
end

function Y = efun(alpha, eqn)
% Y = efun(alpha, eqn)
% This is a function handler that accepts the parameter alpha and returns
% the matrix E_ for solving a second order system re-shaping it as a first
% order system.
Y = [eqn.K_, alpha * eqn.M_; alpha * eqn.M_, eqn.M_];
end
