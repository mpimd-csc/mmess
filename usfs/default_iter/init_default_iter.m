function [result, eqn, opts, oper] = init_default_iter(eqn, opts, oper, flag1, flag2)
% function [result, eqn, opts, oper] = init_default_iter(eqn, opts, oper, flag1, flag2)
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

%% Check input Parameters
% Definition of list of solvers

% Check if eqn.E_ is symmetric and positive definite
if isfield(eqn, 'E_')
    if issymmetric(eqn.E_)
        opts.usfs.default_iter.E_is_sym = true;
        if not(any(diag(eqn.E_ <= 0)))
            opts.usfs.default_iter.E_is_spd = true;
        else
            opts.usfs.default_iter.E_is_spd = false;
        end
    else
        opts.usfs.default_iter.E_is_sym = false;
        opts.usfs.default_iter.E_is_spd = false;
    end
end

% Check if eqn.A_ is symmetric
if isfield(eqn, 'A_')
    if issymmetric(eqn.A_)
        opts.usfs.default_iter.A_is_sym = true;
        if not(any(diag(eqn.A_ <= 0)))
            opts.usfs.default_iter.A_is_spd = true;
        else
            opts.usfs.default_iter.A_is_spd = false;
        end
    else
        opts.usfs.default_iter.A_is_sym = true;
        opts.usfs.default_iter.A_is_spd = false;
    end
end

% Required fields for the iterative solver
if not(isfield(opts, 'usfs')) || not(isfield(opts.usfs, 'default_iter'))

    mess_warn(opts, 'control_data', [' The ''default_iter'' usfs need the ', ...
                                     '''opts.usfs.default_iter'' substructure to be present.']);

end

% Check for the solver used
if isfield(opts.usfs.default_iter, 'method_A')
    if not(exist(opts.usfs.default_iter.method_A, 'file') == 2)
        mess_err(opts, 'control_data', ['Iterative solver method field ''method_A''', ...
                                        ' is an unsupported solver.']);
    end
else
    if opts.usfs.default_iter.A_is_sym
        mess_warn(opts, 'control_data', ['Iterative solver method field ''method_A''', ...
                                         ' is unset. Falling back to PCG.']);
        opts.usfs.default_iter.method_A = 'pcg';
    else
        mess_warn(opts, 'control_data', ['Iterative solver method field ''method_A''', ...
                                         ' is unset. Falling back to GMRES.']);
        opts.usfs.default_iter.method_A = 'gmres';
    end
end

if isfield(opts.usfs.default_iter, 'method_E')
    if not(exist(opts.usfs.default_iter.method_E, 'file') == 2)
        mess_err(opts, 'control_data', ['Iterative solver method field ''method_E''', ...
                                        ' is an unsupported solver.']);
    end
else
    if opts.usfs.default_iter.E_is_spd
        mess_warn(opts, 'control_data', ['Iterative solver method field ''method_E''', ...
                                         ' is unset. Falling back to PCG.']);
        opts.usfs.default_iter.method_E = 'pcg';

    else
        mess_warn(opts, 'control_data', ['Iterative solver method field ''method_E''', ...
                                         ' is unset. Falling back to GMRES.']);
        opts.usfs.default_iter.method_E = 'gmres';
    end
end

% Required residual tolerance for stopping the iterative solver
if isfield(opts.usfs.default_iter, 'res_tol')
    if opts.usfs.default_iter.res_tol < 0
        mess_err(opts, 'control_data', ['Iterative solver residual tolerance value', ...
                                        'is invalid']);
    end
else
    mess_warn(opts, 'control_data', ['Iterative solver residual tolerance value not', ...
                                     ' found. Falling back to default']);
    opts.usfs.default_iter.res_tol = 1e-12;

end

% Number of iterations allowed

if isfield(opts.usfs.default_iter, 'max_iter')
    if opts.usfs.default_iter.max_iter < 0
        mess_err(opts, 'control_data', ['Iterative solver max. iterations', ...
                                        'is invalid']);
    end
else
    mess_warn(opts, 'control_data', ['Iterative solver max. iterations not', ...
                                     ' found. Falling back to default']);
    opts.usfs.default_iter.max_iter = oper.size(eqn, opts, oper);

end

% Restart size for GMRES
if strcmpi(opts.usfs.default_iter.method_A, 'gmres') || ...
        strcmpi(opts.usfs.default_iter.method_E, 'gmres')

    if isfield(opts.usfs.default_iter, 'restIter')
        if opts.usfs.default_iter.restIter < 0
            mess_err(opts, 'control_data', ['GMRES restart iterations value', ...
                                            'is invalid']);
        end
    else
        mess_warn(opts, 'control_data', ['GMRES restart iterations value not', ...
                                         ' found. Falling back to default']);
        opts.usfs.default_iter.restIter = 25;
    end

end

%% Operations
na = nargin;
if isfield(eqn, 'LTV')
    mess_warn(opts, 'not_implemented', ...
              ['''default_iter'' does not yet support', ...
               ' LTV systems.']);
end
if na < 4
    mess_err(opts, 'control_data', 'Number of input Arguments are at least 4');

elseif nargin == 4     % result = init_default_iter(eqn, flag1);
    switch flag1
        case {'A', 'a'}
            [eqn, opts, result] = checkA(eqn, opts);

        case {'E', 'e'}
            [eqn, opts, result] = checkE(eqn, opts);

        otherwise
            mess_err(opts, 'control_data', 'flag1 has to be ''A_'' or ''E_''');
    end

elseif nargin == 5     % result = init_default_iter(eqn,flag1,flag2);
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
            [eqn, result] = checkE(eqn, opts);

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
% This function also defines the preconditioner for A_ by using ICHOL or
% ILU function, depending on whether A_ is symmetric and positive
% definite, or not.

result = isfield(eqn, 'A_');

if result

    result = opts.usfs.default_iter.A_is_sym;

end

if not(isfield(opts.usfs.default_iter, 'PA_R'))

    if not(isfield(opts.usfs.default_iter, 'PA_L'))

        mess_warn(opts, 'control_data', ['No preconditioner for A could be found.', ...
                                         ' Switching to ICHOL/ILU']);
        if opts.usfs.default_iter.A_is_spd
            S = ichol(eqn.A_);
            opts.usfs.default_iter.PA_L = S;
            opts.usfs.default_iter.PA_R = S';

        else
            [L, U] = ilu(-eqn.A_);
            opts.usfs.default_iter.PA_L = L;
            opts.usfs.default_iter.PA_R = U;
        end

    else

        opts.usfs.default_iter.PA_R = [];

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
end

if not(eqn.haveE)
    if isfield(eqn, 'E_')
        mess_err(opts, 'equation_data', ['Detected eqn.E_ where eqn.haveE ' ...
                                         'is 0. You need to set haveE = true or ' ...
                                         'delete E_.']);
    else
        result = true;
    end
else
    result = isfield(eqn, 'E_');
    if result
        result = opts.usfs.default_iter.E_is_spd;
    end

    % Preconditioner for E
    if not(isfield(opts.usfs.default_iter, 'PE_R'))

        if not(isfield(opts.usfs.default_iter, 'PE_L'))

            mess_warn(opts, 'control_data', ['No preconditioner for E could be', ...
                                             ' found. Switching to ICHOL/ILU']);
            if result
                S = ichol(eqn.E_, struct('type', 'nofill', 'michol', 'on'));
                opts.usfs.default_iter.PE_L = S;
                opts.usfs.default_iter.PE_R = S';

            else
                [L, U] = ilu(eqn.E_);
                opts.usfs.default_iter.PE_L = L;
                opts.usfs.default_iter.PE_R = U;
            end
        else

            opts.usfs.default_iter.PE_R = [];

        end

    end

end
end
