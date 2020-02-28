function [result, eqn, opts, oper] = init_dae_2(eqn, opts, oper, flag1, flag2)
% function [result, eqn, opts, oper] = init_dae_2(eqn, opts, oper, flag1, flag2)
% return true or false if Data for A_ and E_ resp. flag1 and flag2  are 
% availabe and correct in eqn.
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
%               2009-2020
%

%% check input Paramters
if(nargin<=3)
    error('MESS:check_data','Number of input Arguments must be at least 3');

%% result = init_dae_2(eqn, opts, oper, flag1);    
elseif(nargin==4)
    switch flag1
        case {'A','a'}
            [eqn,result] = checkA(eqn);
        case {'E','e'}
            [eqn,result] = checkE(eqn);
        otherwise
            error('MESS:check_data','flag1 has to be ''A'' or ''E''');
    end
    
%% result = init_dae_2(eqn, opts, oper,flag1,flag2);
elseif(nargin==5)
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
% A = [ A11 A12;
%       A21   0]
%
if not(isfield(eqn, 'st')) || not(isnumeric(eqn.st))
    error('MESS:st',...
    'Missing or Corrupted st field detected in equation structure.');
end
if not(isfield(eqn,'A_')) || not(isnumeric(eqn.A_))
    error('MESS:equation_data',...
      'Empty or Corrupted field A detected in equation structure.');
end
if  (size(eqn.A_,1) ~= size(eqn.A_,2))
    error('MESS:error_arguments', 'field eqn.A_ has to be quadratic');
end

if(not(issparse(eqn.A_)))
    warning('MESS:check_data','A has to be sparse for best performance');
end
% check if lower right block is empty
if (any(any(eqn.A_(eqn.st+1:end,eqn.st+1:end)))) 
   error('MESS:equation_data',...
      'Corrupted field A detected in equation structure.'); 
end
result = 1;
end

%% checkdata for E_
function [eqn,result] = checkE(eqn)
if not(isfield(eqn, 'haveE')), eqn.haveE = 0; end
if not(isfield(eqn, 'st')) || not(isnumeric(eqn.st))
    error('MESS:st',...
          ['Missing or Corrupted st field detected in equation ' ...
           'structure.']);
end
if eqn.haveE
    if not(isfield(eqn,'E_')) || not(isnumeric(eqn.E_))
        error('MESS:equation_data',...
            'Empty or Corrupted field E detected in equation structure.');
    end
    if  (size(eqn.E_,1) ~= size(eqn.E_,2))
        error('MESS:error_arguments', 'field eqn.E_ has to be quadratic');
    end
    
    if not(issparse(eqn.E_))
        warning('MESS:check_data','E has to be sparse for best performance');
    end
    st = eqn.st;
    % E = [ E1 0;
    %       0  0]
    if full(any([any(eqn.E_(1:st, st + 1:end)), any(eqn.E_(st+1:end,:))]))
        warning('MESS:check_data',['E has to be non-zero only in the ' ...
            'upper left st x st block']);
    end
else
    % E = [ I 0 ]
    %     [ 0 0 ]
    if not(isfield(eqn,'A_')) || not(isnumeric(eqn.A_))
        error('MESS:equation_data',...
            'Empty or Corrupted field A detected in equation structure.');
    end
    st = eqn.st;
    n=size(eqn.A_,1);
    eqn.E_=sparse(1:st,1:st,ones(st, 1),n,n,st);
end
result = 1;
end
