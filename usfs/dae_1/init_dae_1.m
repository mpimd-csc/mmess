function [result, eqn, opts, oper] = init_dae_1(eqn, opts, oper, flag1, flag2)
% function [result, eqn, opts, oper] = init_dae_1(eqn, opts, oper, flag1, flag2)
% return true or false if Data for A_ and E_ resp. flag1 and flag2  are
% available and correct in eqn.
%
%   result = init(eqn,flag1);
%   result = init(eqn,flag1,flag2);
%
%   result = init(eqn,'A')      (==init(eqn,'A','A'));
%   result = init(eqn,'E')      (==init(eqn,'E','E'));
%   result = init(eqn,'A','E')  (==init(eqn,'E','A'));
%
%   Input:
%
%   eqn             structure with data
%   opts            structure containing parameter for the algorithm
%   oper            struct contains usfs for operation with A and E
%   flag1           'A'/'E' to check if A or E is in eqn
%   flag2           'A'/'E' to check if A or E is in eqn
%
%   Output:
%
%   result          1 : data corresponding to flag1 (and flag2) available,
%                   0 : data incomplete
%   eqn             structure with data
%   opts            structure containing parameter for the algorithm
%   oper            struct contains usfs for operation with A and E
%
%   uses no other dae_1 functions
%

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%% check input Parameters
na = nargin;
if not(isfield(eqn, 'LTV'))
    eqn.LTV = false;
end
if na < 3
    mess_err(opts, 'control_data', ...
             'Number of input Arguments are at least 3');

    %% result = init(eqn, flag1);
elseif na == 4
    switch flag1
        case {'A', 'a'}
            if eqn.LTV
                [eqn, result] = checkA_time(eqn, opts);
            else
                [eqn, result] = checkA(eqn, opts);
            end
        case {'E', 'e'}
            if eqn.LTV
                [eqn, result] = checkE_time(eqn, opts);
            else
                [eqn, result] = checkE(eqn, opts);
            end
        otherwise
            mess_err(opts, 'control_data', ...
                     'flag1 has to be ''A_'' or ''E_''');
    end

    %% result = init(eqn,flag1,flag2);
elseif nargin == 5
    switch flag1
        case {'A', 'a'}
            if eqn.LTV
                [eqn, result] = checkA_time(eqn, opts);
            else
                [eqn, result] = checkA(eqn, opts);
            end
            switch flag2
                case {'A', 'a'}
                    if eqn.LTV
                        [eqn, resultA] = checkA_time(eqn, opts);
                    else
                        [eqn, resultA] = checkA(eqn, opts);
                    end
                    result = result && resultA;
                case {'E', 'e'}
                    if eqn.LTV
                        [eqn, resultE] = checkE_time(eqn, opts);
                    else
                        [eqn, resultE] = checkE(eqn, opts);
                    end
                    result = result && resultE;
                otherwise
                    mess_err(opts, 'control_data', ...
                             'flag2 has to be ''A'' or ''E''');
            end
        case {'E', 'e'}
            if eqn.LTV
                [eqn, result] = checkE_time(eqn, opts);
            else
                [eqn, result] = checkE(eqn, opts);
            end
            switch flag2
                case {'A', 'a'}
                    if eqn.LTV
                        [eqn, resultA] = checkA_time(eqn, opts);
                    else
                        [eqn, resultA] = checkA(eqn, opts);
                    end
                    result = result && resultA;
                case {'E', 'e'}
                    if eqn.LTV
                        [eqn, resultE] = checkE_time(eqn, opts);
                    else
                        [eqn, resultE] = checkE(eqn, opts);
                    end
                    result = result && resultE;
                otherwise
                    mess_err(opts, 'control_data', ...
                             'flag2 has to be ''A'' or ''E''');
            end
        otherwise
            mess_err(opts, 'control_data', ...
                     'flag1 has to be ''A'' or ''E''');
    end
end
%% Compute reduced B and C
n = size(eqn.A_, 1);
one = 1:eqn.manifold_dim;
two = eqn.manifold_dim + 1:n;
if not(eqn.LTV)
    if size(eqn.B, 1) > eqn.manifold_dim
        eqn.B = eqn.B(one, :) - eqn.A_(one, two) * ...
                (eqn.A_(two, two) \ eqn.B(two, :));
    end
    if size(eqn.C, 2) > eqn.manifold_dim
        eqn.C = eqn.C(:, one) - ...
                (eqn.C(:, two) / eqn.A_(two, two)) * eqn.A_(two, one);
    end
end
end

%% checkdata for A_
function [eqn, result] = checkA(eqn, opts)

if not(isfield(eqn, 'A_')) || not(isnumeric(eqn.A_))
    mess_err(opts, 'equation_data', ...
             'Empty or Corrupted field A detected in equation structure.');
end

if not(size(eqn.A_, 1) == size(eqn.A_, 2))
    mess_err(opts, 'error_arguments', 'field eqn.A_ has to be quadratic');
end

if not(issparse(eqn.A_))
    mess_warn(opts, 'control_data', 'A is not sparse');
end

if not(isfield(eqn, 'manifold_dim'))    || not(isnumeric(eqn.manifold_dim))
    mess_err(opts, 'equations_data', ...
             ['Missing or corrupted manifold_dim field detected in ' ...
              'equation structure.']);
end
result = true;
end

%% checkdata for E_
function [eqn, result] = checkE(eqn, opts)
if not(isfield(eqn, 'haveE'))
    eqn.haveE = false;
end
if not(isfield(eqn, 'manifold_dim')) || not(isnumeric(eqn.manifold_dim))
    mess_err(opts, 'equation_data', ...
             ['Missing or corrupted manifold_dim field detected in ' ...
              'equation structure.']);
elseif (not(isfield(eqn, 'E_')) || not(isnumeric(eqn.E_))) && eqn.haveE
    mess_err(opts, 'equation_data', ...
             'Empty or Corrupted field E detected in equation structure.');
end
n_ode = eqn.manifold_dim;
if not(isfield(eqn, 'A_')) || not(isnumeric(eqn.A_))
    mess_err(opts, 'equation_data', ...
             'Empty or Corrupted field A detected in equation structure.');
end
n = size(eqn.A_, 1);

if not(eqn.haveE)
    if isfield(eqn, 'E_')
        mess_err(opts, 'equation_data', ...
                 ['Detected eqn.E_ where eqn.haveE ' ...
                  'is 0. You need to set haveE = true or delete E_.']);
    else
        result = true;
    end
else
    if not(size(eqn.E_, 1) == size(eqn.E_, 2))
        mess_err(opts, 'error_arguments', ...
                 'field eqn.E_ has to be quadratic');
    end
    if not(issparse(eqn.E_))
        mess_warn(opts, 'control_data', ...
                  'E is not sparse');
    end
    % check size(A) == size(E)?
    if not(n == size(eqn.E_, 1))
        mess_err(opts, 'error_arguments', ...
                 'dimensions of E and A must coincide');
    end
    % E = [ E1 0 ]
    %     [ 0  0 ]
    if full(any([any(eqn.E_(1:n_ode, n_ode + 1:end)), any(eqn.E_(n_ode + 1:end, :))]))
        mess_warn(opts, 'control_data', ...
                  'E has to be non-zero only in n_ode x n_ode block');
    end
    % result: bool; without 'full()' result: 1x1 sparse
    result = true;
end
end

%% checkdata for A_
function [eqn, result] = checkA_time(eqn, opts)
if not(isfield(eqn, 'A_time')) || not(isa(eqn.A_time, 'function_handle'))
    mess_err(opts, 'equation_data', ...
             'Empty or Corrupted field A_time detected in equation structure.');
end
A = eqn.A_time(opts.t0);
if not(isnumeric(A))
    mess_err(opts, 'equation_data', ...
             ['Empty or Corrupted field eqn.A_time(t) detected in', ...
              ' equation structure.']);
end
if not(size(A, 1) == size(A, 2))
    mess_err(opts, 'error_arguments', ...
             'field eqn.A_time(t) has to be quadratic');
end
if not(issparse(A))
    mess_warn(opts, 'control_data', 'A is not sparse');
end
if not(isfield(eqn, 'manifold_dim')) || not(isnumeric(eqn.manifold_dim))
    mess_err(opts, 'equations_data', ...
             ['Missing or corrupted manifold_dim field detected in ' ...
              'equation structure.']);
end
result = true;
end

%% checkdata for E_
function [eqn, result] = checkE_time(eqn, opts)
if not(isfield(eqn, 'haveE'))
    eqn.haveE = false;
end
if not(isfield(eqn, 'manifold_dim')) || not(isnumeric(eqn.manifold_dim))
    mess_err(opts, 'equation_data', ...
             ['Missing or corrupted manifold_dim field detected in ' ...
              'equation structure.']);
end
n_ode = eqn.manifold_dim;
if not(isfield(eqn, 'A_time')) || not(isa(eqn.A_time, 'function_handle'))
    mess_err(opts, 'equation_data', ...
             'Empty or Corrupted field A_time detected in equation structure.');
end
A = eqn.A_time(opts.t0);
if not(isnumeric(A))
    mess_err(opts, 'equation_data', ...
             ['Empty or Corrupted field A_time(t) detected in ', ...
              'equation structure.']);
end
n = size(A, 1);

if not(eqn.haveE)
    if isfield(eqn, 'E_time')
        mess_err(opts, 'equation_data', ...
                 ['Detected eqn.E_time where eqn.haveE '...
                  'is 0. You need to set haveE = true or delete E_']);
    else
        result = true;
    end
else
    if not(isfield(eqn, 'E_time')) || not(isa(eqn.E_time, 'function_handle'))
        mess_err(opts, 'equation_data', ...
                 ['Empty or Corrupted field E_time detected in ', ...
                  'equation structure.']);
    end
    E = eqn.E_time(opts.t0);
    if not(isnumeric(E))
        mess_err(opts, 'equation_data', ...
                 ['Empty or Corrupted field eqn.E_time(t) detected in ', ...
                  'equation structure.']);
    end
    if not(size(E, 1) == size(E, 2))
        mess_err(opts, 'error_arguments', ...
                 'field eqn.E_time(t) has to be quadratic');
    end
    if not(issparse(E))
        mess_warn(opts, 'control_data', ...
                  'eqn.E_time(t) is not sparse');
    end
    % check size(A) == size(E)?
    if not(n == size(E, 1))
        mess_err(opts, 'error_arguments', ...
                 'dimensions of eqn.E_time(t) and eqn.A_time(t) must coincide');
    end
    % E = [ E1 0 ]
    %     [ 0  0 ]
    if full(any([any(E(1:n_ode, n_ode + 1:end)), any(E(n_ode + 1:end, :))]))
        mess_warn(opts, 'control_data', ...
                  'eqn.E_time(t) has to be non-zero only in n_ode x n_ode block');
    end
    % result: bool; without 'full()' result: 1x1 sparse
    result = true;
end
end
