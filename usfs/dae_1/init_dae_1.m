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
% This file is part of the M-M.E.S.S. project 
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%


%% check input Parameters
na = nargin;
if not(isfield(eqn, 'LTV')),  eqn.LTV=0; end
if(na<3)
    error('MESS:control_data','Number of input Arguments are at least 3');

%% result = init(eqn, flag1);
elseif(na==4)
    switch flag1
        case {'A','a'}
            if eqn.LTV
                [eqn,result] = checkA_time(eqn,opts);
            else
                [eqn,result] = checkA(eqn);
            end
        case {'E','e'}
            if eqn.LTV
                [eqn,result] = checkE_time(eqn,opts);
            else
                [eqn,result] = checkE(eqn);
            end
        otherwise
            error('MESS:control_data','flag1 has to be ''A_'' or ''E_''');
    end

%% result = init(eqn,flag1,flag2);
elseif(nargin==5)
    switch flag1
        case {'A','a'}
            if eqn.LTV
                [eqn, result] = checkA_time(eqn,opts);
            else
                [eqn, result] = checkA(eqn);
            end
            switch flag2
                case {'A','a'}
                    if eqn.LTV
                        [eqn, resultA] = checkA_time(eqn,opts);
                    else
                        [eqn, resultA] = checkA(eqn);
                    end
                    result = result && resultA;
                case {'E','e'}
                    if eqn.LTV
                        [eqn, resultE] = checkE_time(eqn,opts);
                    else
                        [eqn, resultE]= checkE(eqn);
                    end
                    result = result && resultE;
                otherwise
                    error('MESS:control_data', ...
                          'flag2 has to be ''A'' or ''E''');
            end
        case {'E','e'}
            if eqn.LTV
                [eqn, result] = checkE_time(eqn,opts);
            else
                [eqn, result] = checkE(eqn);
            end
            switch flag2
                case {'A','a'}
                    if eqn.LTV
                        [eqn, resultA] = checkA_time(eqn,opts);
                    else
                        [eqn, resultA] = checkA(eqn);
                    end
                    result = result && resultA;
                case {'E','e'}
                    if eqn.LTV
                        [eqn, resultE] = checkE_time(eqn,opts);
                    else
                        [eqn, resultE]= checkE(eqn);
                    end
                    result = result && resultE;
                otherwise
                    error('MESS:control_data', ...
                          'flag2 has to be ''A'' or ''E''');
            end
        otherwise
            error('MESS:control_data',...
                  'flag1 has to be ''A'' or ''E''');
    end
end
%% Compute reduced B and C
n = size(eqn.A_,1);
st = eqn.st;
one = 1:st;
two = st + 1 : n;
if not(eqn.LTV)
    if size(eqn.B, 1) > st
        eqn.B = eqn.B(one, :) - eqn.A_(one, two) ...
            * (eqn.A_(two, two) \ eqn.B(two, :));
    end
    if size(eqn.C, 2) > st
        eqn.C = eqn.C( : , one) - (eqn.C( : , two) ...
            / eqn.A_(two, two)) * eqn.A_(two, one);
    end
end
end

%% checkdata for A_
function [eqn,result] = checkA(eqn)
if not(isfield(eqn,'A_')) || not(isnumeric(eqn.A_))
    error('MESS:equation_data',...
        'Empty or Corrupted field A detected in equation structure.');
end
if  (size(eqn.A_,1) ~= size(eqn.A_,2))
    error('MESS:error_arguments', 'field eqn.A_ has to be quadratic');
end
if(not(issparse(eqn.A_)))
    warning('MESS:control_data','A is not sparse');
end
if not(isfield(eqn, 'st'))    || not(isnumeric(eqn.st))
    error('MESS:st',...
        'Missing or Corrupted st field detected in equation structure.');
end
result = 1;
end

%% checkdata for E_
function [eqn,result] = checkE(eqn)
if not(isfield(eqn, 'haveE')), eqn.haveE = 0; end
if not(isfield(eqn, 'st'))    || not(isnumeric(eqn.st))
    error('MESS:st',...
        'Missing or Corrupted st field detected in equation structure.');
elseif (not(isfield(eqn,'E_')) || not(isnumeric(eqn.E_))) && eqn.haveE
    error('MESS:equation_data',...
        'Empty or Corrupted field E detected in equation structure.');
end
st = eqn.st;
if not(isfield(eqn,'A_')) || not(isnumeric(eqn.A_))
    error('MESS:equation_data',...
        'Empty or Corrupted field A detected in equation structure.');
end
n=size(eqn.A_,1);

if not(eqn.haveE)
    if isfield(eqn, 'E_')
        error('MESS:equation_data', ['Detected eqn.E_ where eqn.haveE ' ...
            'is 0. You need to set haveE=1 or delete E_.']);
    else
        result = 1;
    end
else
    if  (size(eqn.E_,1) ~= size(eqn.E_,2))
        error('MESS:error_arguments', 'field eqn.E_ has to be quadratic');
    end
    if(not(issparse(eqn.E_)))
        warning('MESS:control_data','E is not sparse');
    end
    % check size(A) == size(E)?
    if (n~=size(eqn.E_,1))
        error('MESS:error_arguments',...
            'dimensions of E and A must coincide');
    end
    % E = [ E1 0 ]
    %     [ 0  0 ]
    if full(any([any(eqn.E_(1:st, st + 1:end)), any(eqn.E_(st+1:end,:))]))
        warning('MESS:control_data',...
            'E has to be non-zero only in st x st block');
    end
    % result: bool; without 'full()' result: 1x1 sparse
    result = 1;
end
end

%% checkdata for A_
function [eqn, result] = checkA_time(eqn,opts)
if not(isfield(eqn, 'A_time')) || not(isa(eqn.A_time,'function_handle'))
    error('MESS:equation_data',...
        'Empty or Corrupted field A_time detected in equation structure.');
end
A = eqn.A_time(opts.t0);
if not(isnumeric(A))
    error('MESS:equation_data',...
        'Empty or Corrupted field eqn.A_time(t) detected in equation structure.');
end
if  (size(A,1) ~= size(A,2))
    error('MESS:error_arguments', 'field eqn.A_time(t) has to be quadratic');
end
if(not(issparse(A)))
    warning('MESS:control_data','A is not sparse');
end
if not(isfield(eqn, 'st')) || not(isnumeric(eqn.st))
    error('MESS:st',...
        'Missing or Corrupted st field detected in equation structure.');
end
result = 1;
end

%% checkdata for E_
function [eqn, result] = checkE_time(eqn,opts)
if not(isfield(eqn, 'haveE')), eqn.haveE = 0; end
if not(isfield(eqn, 'st')) || not(isnumeric(eqn.st))
    error('MESS:st',...
        'Missing or Corrupted st field detected in equation structure.');
end
st = eqn.st;
if not(isfield(eqn, 'A_time')) || not(isa(eqn.A_time,'function_handle'))
    error('MESS:equation_data',...
        'Empty or Corrupted field A_time detected in equation structure.');
end
A = eqn.A_time(opts.t0);
if not(isnumeric(A))
    error('MESS:equation_data',...
        'Empty or Corrupted field A_time(t) detected in equation structure.');
end
n=size(A,1);

if not(eqn.haveE)
    if isfield(eqn, 'E_time')
        error('MESS:equation_data',['Detected eqn.E_time where eqn.haveE '...
            'is 0. You need to set haveE=1 or delete E_']);
    else
        result = 1;
    end
else
    if not(isfield(eqn, 'E_time')) || not(isa(eqn.E_time,'function_handle'))
        error('MESS:equation_data',...
            'Empty or Corrupted field E_time detected in equation structure.');
    end
    E = eqn.E_time(opts.t0);
    if not(isnumeric(E))
        error('MESS:equation_data',...
            'Empty or Corrupted field eqn.E_time(t) detected in equation structure.');
    end
    if  (size(E,1) ~= size(E,2))
        error('MESS:error_arguments', 'field eqn.E_time(t) has to be quadratic');
    end
    if(not(issparse(E)))
        warning('MESS:control_data','eqn.E_time(t) is not sparse');
    end
    % check size(A) == size(E)?
    if (n~=size(E,1))
        error('MESS:error_arguments',...
            'dimensions of eqn.E_time(t) and eqn.A_time(t) must coincide');
    end
    % E = [ E1 0 ]
    %     [ 0  0 ]
    if full(any([any(E(1:st, st + 1:end)), any(E(st+1:end,:))]))
        warning('MESS:control_data',...
            'eqn.E_time(t) has to be non-zero only in st x st block');
    end
    % result: bool; without 'full()' result: 1x1 sparse
    result = 1;
end
end
