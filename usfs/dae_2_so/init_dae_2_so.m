function [result, eqn, opts, oper] = init_dae_2_so(eqn, opts, oper, flag1, flag2)
% function [result, eqn, opts, oper] = init_dae_2_so(eqn, opts, oper, flag1, flag2)
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
%   result             1 if data corresponding to flag1 (and flag2) are available , 0 data are not available
%   eqn             structure with data
%   opts            structure containing parameter for the algorithm
%   oper            struct contains function handles for operation with A and E

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%



%% check input Parameters
na = nargin;
if(na<=3)
    error('MESS:check_data','Number of input Arguments must be at least 3');

    %% result = init_dae_2(eqn, flag1);
elseif(na==4)
    switch flag1
        case {'A','a'}
            [eqn,result] = checkA(eqn);
        case {'E','e'}
            [eqn,result] = checkE(eqn);
        otherwise
            error('MESS:check_data','flag1 has to be ''A'' or ''E''');
    end

    %% result = init_dae_2(eqn,flag1,flag2);
elseif(na==5)
    switch flag1
        case {'A','a'}
            [eqn,result] = checkA(eqn);
            switch flag2
                case {'A','a'}
                    [eqn,resultA] = checkA(eqn);
                    result = result && resultA;
                case {'E','e'}
                    [eqn,resultE] = checkE(eqn);
                    result = result &&resultE;
                otherwise
                    error('MESS:check_data','flag2 has to be ''A'' or ''E''');
            end
        case {'E','e'}
            [eqn, result] = checkE(eqn);
            switch flag2
                case {'A','a'}
                    [eqn,resultA] = checkA(eqn);
                    result = result && resultA;
                case {'E','e'}
                    [eqn,resultE] = checkE(eqn);
                    result =  result && resultE;
                otherwise
                    error('MESS:check_data','flag2 has to be ''A'' or ''E''');
            end
        otherwise
            error('MESS:check_data','flag1 has to be ''A'' or ''E''');
    end
end
end

%% checkdata for A_
function [eqn,result] = checkA(eqn)

if not(isfield(eqn,'E_')) || not(isnumeric(eqn.E_))
    error('MESS:equation_data',...
        'Empty or Corrupted field D detected in equation structure.');
elseif not(issparse(eqn.E_))
    warning('MESS:control_data','D is not sparse');
end

if not(isfield(eqn,'K_')) || not(isnumeric(eqn.K_))
    error('MESS:equation_data',...
        'Empty or Corrupted field K detected in equation structure.');
elseif not(issparse(eqn.K_))
    warning('MESS:control_data','K is not sparse');
end

if not(isfield(eqn,'G_')) || not(isnumeric(eqn.G_))
    error('MESS:equation_data',...
        'Empty or Corrupted field K detected in equation structure.');
elseif not(issparse(eqn.G_))
    warning('MESS:control_data','G is not sparse');
end

if  (size(eqn.E_,1) ~= size(eqn.E_,2))
    error('MESS:error_arguments', 'field eqn.E_ has to be quadratic');
end

if  (size(eqn.K_,1) ~= size(eqn.K_,2))
    error('MESS:error_arguments', 'field eqn.K_ has to be quadratic');
end

if  (size(eqn.E_,1) ~= size(eqn.G_,2))
    error('MESS:error_arguments', 'field eqn.G_ has invalid number of columns');
end

result = 1;
end

%% checkdata for E_
function [eqn,result] = checkE(eqn)
if not(isfield(eqn, 'haveE')), eqn.haveE = 1; end

if not(isfield(eqn,'M_')) || not(isnumeric(eqn.M_))
    error('MESS:equation_data',...
        'Empty or Corrupted field M detected in equation structure.');
end
if  (size(eqn.M_,1) ~= size(eqn.M_,2))
    error('MESS:error_arguments', 'field eqn.M_ has to be quadratic');
end
if(not(issparse(eqn.M_)))
    warning('MESS:check_data','M is not sparse');
end

if not(isfield(eqn, 'alpha')) || not(isnumeric(eqn.alpha))
   error('MESS:equation_data',...
         'No parameter alpha given for shifting infinite eigenvalues of the pencil');
end
result=1;
end
