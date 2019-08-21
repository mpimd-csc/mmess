function [result, eqn, opts, oper] = init_dae_1_so(eqn, opts, oper, flag1, flag2)
% function [result, eqn, opts, oper] = init_dae_1_so(eqn, opts, oper, flag1, flag2)
% return true or false if Data for A and E resp. flag1 and flag2  are availabe
% and correct in eqn.
%
%   result = init_so_1(eqn,flag1);
%   result = init_so_1(eqn,flag1,flag2);
%
%   result = init_so_1(eqn,'A')    (==init_so_1(eqn,'A','A'));
%   result = init_so_1(eqn,'E')    (==init_so_1(eqn,'E','E'));
%   result = init_so_1(eqn,'A','E')  (==init_so_1(eqn,'E','A'));
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
%   uses no other dae_1_so function

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

%% start checking
na = nargin;
if(na<=3)
    error('MESS:control_data','Number of input Arguments are at least 3');

%% result = init_so_1(eqn, flag1);    
elseif(na==4)
    switch flag1
        case {'A','a'}
            [eqn,result] = checkA(eqn);
        case {'E','e'}
            [eqn,result] = checkE(eqn);
        otherwise
            error('MESS:control_data','flag1 has to be ''A'' or ''E''');
    end
    
%% result = init_so_1(eqn,flag1,flag2);
elseif(na==5)
    switch flag1
        case {'A','a'}
            [eqn,result] = checkA(eqn);
            switch flag2
                case {'A','a'}
                    [eqn,resultA] = checkA(eqn);
                    result = result && resultA;
                case {'E','e'}
                    [eqn, resultE] = checkE(eqn);
                    result =  result && resultE;
                otherwise
                    error('MESS:control_data','flag2 has to be ''A'' or ''E''');
            end
        case {'E','e'}
            [eqn,result] = checkE(eqn);
            switch flag2
                case {'A','a'}
                    [eqn,resultA] = checkA(eqn);
                    result =  result && resultA;
                case {'E','e'}
                    [eqn,resultE] = checkE(eqn);
                    result = result && resultE;
                otherwise
                    error('MESS:control_data','flag2 has to be ''A'' or ''E''');
            end
        otherwise
            error('MESS:control_data','flag1 has to be ''A'' or ''E''');
    end 
end
end

%% checkdata for A

function [eqn,result] = checkA(eqn)
if (not(isfield(eqn,'K_')) || not(isnumeric(eqn.K_)))
    error('MESS:equation_data',...
        'Empty or Corrupted field K detected in equation structure.')
end
if(not(issparse(eqn.K_)))
    warning('MESS:control_data','K is not sparse');
end
[n1k, n2k] = size(eqn.K_);
if n1k ~= n2k
    error('MESS:equation_data',...
        'K has to be quadratic')
end
result = 1;

end

%% checkdata for E
function [eqn,result] = checkE(eqn)
if (not(isfield(eqn,'M_')) || not(isnumeric(eqn.M_)))
    error('MESS:equation_data',...
        'Empty or Corrupted field M detected in equation structure.')
elseif (not(isfield(eqn,'E_')) || not(isnumeric(eqn.E_)))
    error('MESS:equation_data',...
        'Empty or Corrupted field D detected in equation structure.')
end
if not(isfield(eqn, 'nd'))    || not(isnumeric(eqn.nd))
    error('MESS:nd',...
    'Missing or Corrupted nd field detected in equation structure.');
end
if(not(issparse(eqn.M_)))
    warning('MESS:control_data','M is not sparse');
end
if(not(issparse(eqn.E_)))
    warning('MESS:control_data','D is not sparse');
end
nd = eqn.nd; 
[n1m, n2m] = size(eqn.M_);
[n1d, n2d] = size(eqn.E_);
if n1m ~= n2m
    error('MESS:equation_data',...
        'M has to be quadratic')
end
if n1d ~= n2d
    error('MESS:equation_data',...
        'D has to be quadratic')
end
if n1m ~= n1d
    error('MESS:equation_data',...
        'M and D must have same size')
end
if full(any([any(eqn.M_(1:nd, nd + 1:end)), any(eqn.M_(nd+1:end,:))]))
    warning('MESS:control_data','M has to be non-zero only in nd x nd block');
end
if full(any([any(eqn.E_(1:nd, nd + 1:end)), any(eqn.E_(nd+1:end,:))]))
    warning('MESS:control_data','D has to be non-zero only in nd x nd block');
end
result = 1;

end
