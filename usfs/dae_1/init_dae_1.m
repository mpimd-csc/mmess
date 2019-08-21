function [result, eqn, opts, oper] = init_dae_1(eqn, opts, oper, flag1, flag2)
% function [result, eqn, opts, oper] = init_dae_1(eqn, opts, oper, flag1, flag2)
% return true or false if Data for A_ and E_ resp. flag1 and flag2  are a
% vailabe and correct in eqn.
%
%   result = init(eqn,flag1);
%   result = init(eqn,flag1,flag2);
%
%   result = init(eqn,'A')    (==init(eqn,'A','A'));
%   result = init(eqn,'E')    (==init(eqn,'E','E'));
%   result = init(eqn,'A','E')  (==init(eqn,'E','A'));
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
%   uses no other dae_1 functions

%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, see <http://www.gnu.org/licenses/>.
%
% Copyright (C) Jens Saak, Martin Koehler, Peter Benner and others
%               2009-2019
%

%% check input Paramters
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
                    error('MESS:control_data','flag2 has to be ''A'' or ''E''');
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
                    error('MESS:control_data','flag2 has to be ''A'' or ''E''');
            end
        otherwise
            error('MESS:control_data','flag1 has to be ''A'' or ''E''');
    end
end
%% Compute reduced B and C
st = eqn.st;
if not(eqn.LTV)
    if size(eqn.B, 1) > st
        eqn.B = eqn.B(1 : st, :) - eqn.A_(1 : st, st + 1 : end) ...
            * (eqn.A_(st + 1 : end, st + 1 : end) \ eqn.B(st + 1 : end, :));
    end
    if size(eqn.C, 2) > st
        eqn.C = eqn.C( : , 1 : st) - (eqn.C( : , st + 1 : end) ...
            / eqn.A_(st +1 : end, st + 1 : end)') * eqn.A_(1 : st, st+1 : end)';
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
    % E = [ I 0 ]
    %     [ 0 0 ]
    eqn.E_=sparse(1:st,1:st,ones(st, 1),n,n,st);
else
    if  (size(eqn.E_,1) ~= size(eqn.E_,2))
        error('MESS:error_arguments', 'field eqn.E_ has to be quadratic');
    end
    if(not(issparse(eqn.E_)))
        warning('MESS:control_data','E is not sparse');
    end
    % check size(A) == size(E)?
    if (n~=size(eqn.E_,1))
        error('MESS:error_arguments','dimensions of E and A must coincide');
    end
    % E = [ E1 0 ]
    %     [ 0  0 ]
    if full(any([any(eqn.E_(1:st, st + 1:end)), any(eqn.E_(st+1:end,:))]))
        warning('MESS:control_data','E has to be non-zero only in st x st block');
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
    % E = [ I 0 ]
    %     [ 0 0 ]
    eqn.E_=sparse(1:st,1:st,ones(st, 1),n,n,st);
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
