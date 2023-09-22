function [result, eqn, opts, oper] = init_dae_2(eqn, opts, oper, flag1, flag2)
% function [result, eqn, opts, oper] = init_dae_2(eqn, opts, oper, flag1, flag2)
% return true or false if Data for A_ and E_ resp. flag1 and flag2  are
% available and correct in eqn.
%
%   result = init_dae_2(eqn,flag1);
%   result = init_dae_2(eqn,flag1,flag2);
%
%   result = init_dae_2(eqn,'A')    (==init_dae_2(eqn,'A','A'));
%   result = init_dae_2(eqn,'E')    (==init_dae_2(eqn,'E','E'));
%   result = init_dae_2(eqn,'A','E')  (==init_dae_2(eqn,'E','A'));
%
%   Input:
%
%   eqn             structure with data
%   opts            structure containing parameter for the algorithm
%   oper            struct contains function handles for operation with A and E
%   flag1           'A'/'E' to check if A or E is in eqn
%   flag2           'A'/'E' to check if A or E is in eqn
%
%   Output:
%
%   result             1 if data corresponding to flag1 (and flag2)
%                      are available ,
%                   0 data are not available
%   eqn             structure with data
%   opts            structure containing parameter for the algorithm
%   oper            struct contains function handles for operation with A and E
%
%   uses no other dae_2 functions

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%% check input Parameters
if nargin <= 3
    mess_err(opts, 'check_data', 'Number of input Arguments must be at least 3');

    %% result = init_dae_2(eqn, opts, oper, flag1);
elseif nargin == 4
    switch flag1
        case {'A', 'a'}
            [eqn, result] = checkA(eqn, opts);
        case {'E', 'e'}
            [eqn, result] = checkE(eqn, opts);
        otherwise
            mess_err(opts, 'check_data', 'flag1 has to be ''A'' or ''E''');
    end

    %% result = init_dae_2(eqn, opts, oper,flag1,flag2);
elseif nargin == 5
    switch flag1
        case {'A', 'a'}
            [eqn, result] = checkA(eqn, opts);
            switch flag2
                case {'A', 'a'}
                    [eqn, resultA] = checkA(eqn, opts);
                    result = result && resultA;
                case {'E', 'e'}
                    [eqn, resultE] = checkE(eqn, opts);
                    result = result && resultE;
                otherwise
                    mess_err(opts, 'check_data', ...
                             'flag2 has to be ''A'' or ''E''');
            end
        case {'E', 'e'}
            [eqn, result] = checkE(eqn, opts);
            switch flag2
                case {'A', 'a'}
                    [eqn, resultA] = checkA(eqn, opts);
                    result = result && resultA;
                case {'E', 'e'}
                    [eqn, resultE] = checkE(eqn, opts);
                    result =  result && resultE;
                otherwise
                    mess_err(opts, 'check_data', ...
                             'flag2 has to be ''A'' or ''E''');
            end
        otherwise
            mess_err(opts, 'check_data', 'flag1 has to be ''A'' or ''E''');
    end
end

end

%% checkdata for A_
function [eqn, result] = checkA(eqn, opts)
% A = [ A11 A12;
%       A21   0]
%
if not(isfield(eqn, 'manifold_dim')) || not(isnumeric(eqn.manifold_dim))
    mess_err(opts, 'error_arguments', ...
             ['Missing or corrupted manifold_dim field detected in ' ...
              'equation structure.']);
end
if not(isfield(eqn, 'A_')) || not(isnumeric(eqn.A_))
    mess_err(opts, 'equation_data', ...
             'Empty or corrupted field A detected in equation structure.');
end
if not(size(eqn.A_, 1) == size(eqn.A_, 2))
    mess_err(opts, 'error_arguments', 'field eqn.A_ has to be quadratic');
end

if not(issparse(eqn.A_))
    mess_warn(opts, 'check_data', 'A has to be sparse for best performance');
end
% check if lower right block is empty
if any(any(eqn.A_(eqn.manifold_dim + 1:end, eqn.manifold_dim + 1:end)))
    mess_err(opts, 'equation_data', ...
             'Corrupted field A detected in equation structure.');
end
result = true;
end

%% checkdata for E_
function [eqn, result] = checkE(eqn, opts)
if not(isfield(eqn, 'haveE'))
    eqn.haveE = false;
end
if not(isfield(eqn, 'manifold_dim')) || not(isnumeric(eqn.manifold_dim))
    mess_err(opts, 'error_arguments', ...
             ['Missing or corrupted manifold_dim field detected in equation ' ...
              'structure.']);
end
if eqn.haveE
    if not(isfield(eqn, 'E_')) || not(isnumeric(eqn.E_))
        mess_err(opts, 'equation_data', ...
                 'Empty or corrupted field E detected in equation structure.');
    end
    if not(size(eqn.E_, 1) == size(eqn.E_, 2))
        mess_err(opts, 'error_arguments', 'field eqn.E_ has to be quadratic');
    end

    if not(issparse(eqn.E_))
        mess_warn(opts, 'check_data', ...
                  'E has to be sparse for best performance');
    end
    n_ode = eqn.manifold_dim;
    % E = [ E1 0;
    %       0  0]
    if full(any([any(eqn.E_(1:n_ode, n_ode + 1:end)), any(eqn.E_(n_ode + 1:end, :))]))
        mess_warn(opts, 'check_data', ['E has to be non-zero only in the ' ...
                                       'upper left n_ode x n_ode block']);
    end
else
    % E = [ I 0 ]
    %     [ 0 0 ]
    if not(isfield(eqn, 'A_')) || not(isnumeric(eqn.A_))
        mess_err(opts, 'equation_data', ...
                 'Empty or corrupted field A detected in equation structure.');
    end
    n_ode = eqn.manifold_dim;
    n = size(eqn.A_, 1);
    eqn.E_ = sparse(1:n_ode, 1:n_ode, ones(n_ode, 1), n, n, n_ode);
end
result = true;
end
