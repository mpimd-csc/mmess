function [result, eqn, opts, oper] = init_dae_1_so(eqn, opts, oper, flag1, flag2)
% function [result, eqn, opts, oper] = init_dae_1_so(eqn, opts, oper, flag1, flag2)
% return true or false if data for A and E resp. flag1 and flag2  are available
% and correct in eqn.
%
%   result = init_dae_1_so(eqn,flag1);
%   result = init_dae_1_so(eqn,flag1,flag2);
%
%   result = init_dae_1_so(eqn,'A')    (==init_dae_1_so(eqn,'A','A'));
%   result = init_dae_1_so(eqn,'E')    (==init_dae_1_so(eqn,'E','E'));
%   result = init_dae_1_so(eqn,'A','E')  (==init_dae_1_so(eqn,'E','A'));
%
%   Input:
%
%   eqn             structure with data
%   opts            structure containing parameter for the algorithm
%   oper            struct contains function handles for operation with A and E
%   flag1           'A' or 'E' to check if A or E is represented
%                   correctly in eqn
%   flag2           'A' or 'E' to check if A or E is represented
%                   correctly in eqn
%
%   Output:
%
%   result          1 if data corresponding to flag1 (and flag2)
%                   are available , 0 data are not available
%   eqn             structure with data
%   opts            structure containing parameter for the algorithm
%   oper            struct contains function handles for operation with A and E
%
%   uses no other dae_1_so function

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%% start checking
if nargin <= 3
    mess_err(opts, 'control_data', 'Number of input Arguments are at least 3');

    %% result = init_dae_1_so_1(eqn, flag1);
elseif nargin == 4
    switch flag1
        case {'A', 'a'}
            [eqn, result] = checkA(eqn);
        case {'E', 'e'}
            [eqn, result] = checkE(eqn);
        otherwise
            mess_err(opts, 'control_data', 'flag1 has to be ''A'' or ''E''');
    end

    %% result = init_dae_1_so_1(eqn,flag1,flag2);
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
                    result =  result && resultE;
                otherwise
                    mess_err(opts, 'control_data', ...
                             'flag2 has to be ''A'' or ''E''');
            end
        case {'E', 'e'}
            [eqn, result] = checkE(eqn, opts);
            switch flag2
                case {'A', 'a'}
                    [eqn, resultA] = checkA(eqn, opts);
                    result =  result && resultA;
                case {'E', 'e'}
                    [eqn, resultE] = checkE(eqn, opts);
                    result = result && resultE;
                otherwise
                    mess_err(opts, 'control_data', ...
                             'flag2 has to be ''A'' or ''E''');
            end
        otherwise
            mess_err(opts, 'control_data', 'flag1 has to be ''A'' or ''E''');
    end
end

if not(isfield(eqn, 'haveE'))
    eqn.haveE = true;
end

end

%% checkdata for A

function [eqn, result] = checkA(eqn, opts)
if not(isfield(eqn, 'K_')) || not(isnumeric(eqn.K_))
    mess_err(opts, 'equation_data', ...
             'Empty or Corrupted field K detected in equation structure.');
end
if not(issparse(eqn.K_))
    mess_warn(opts, 'control_data', 'K is not sparse');
end
[n1k, n2k] = size(eqn.K_);
if not(n1k == n2k)
    mess_err(opts, 'equation_data', ...
             'K has to be quadratic');
end
result = true;

end

%% checkdata for E
function [eqn, result] = checkE(eqn, opts)
if not(isfield(eqn, 'M_')) || not(isnumeric(eqn.M_))
    mess_err(opts, 'equation_data', ...
             'Empty or Corrupted field M detected in equation structure.');
elseif not(isfield(eqn, 'E_')) || not(isnumeric(eqn.E_))
    mess_err(opts, 'equation_data', ...
             'Empty or Corrupted field D detected in equation structure.');
end
if not(isfield(eqn, 'manifold_dim')) || not(isnumeric(eqn.manifold_dim))
    mess_err(opts, 'equation_data', ...
             ['Missing or corrupted manifold_dim field detected in ' ...
              'equation structure.']);
end
if not(issparse(eqn.M_))
    mess_warn(opts, 'control_data', 'M is not sparse');
end
if not(issparse(eqn.E_))
    mess_warn(opts, 'control_data', 'D is not sparse');
end
manifold_dim = eqn.manifold_dim;
[n1m, n2m] = size(eqn.M_);
[n1d, n2d] = size(eqn.E_);
if not(n1m == n2m)
    mess_err(opts, 'equation_data', ...
             'M has to be quadratic');
end
if not(n1d == n2d)
    mess_err(opts, 'equation_data', ...
             'D has to be quadratic');
end
if not(n1m == n1d)
    mess_err(opts, 'equation_data', ...
             'M and D must have same size');
end
if full(any([any(eqn.M_(1:manifold_dim, manifold_dim + 1:manifold_dim)), ...
             any(eqn.M_(manifold_dim + 1:end, :))]))
    mess_warn(opts, 'control_data', ...
              'M has to be non-zero only in manifold_dim x manifold_dim block');
end
if full(any([any(eqn.E_(1:manifold_dim, manifold_dim + 1:end)), ...
             any(eqn.E_(manifold_dim + 1:end, :))]))
    mess_warn(opts, 'control_data', ...
              'D has to be non-zero only in manifold_dim x manifold_dim block');
end
result = true;

end
