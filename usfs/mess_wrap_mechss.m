function [eqn, oper] = mess_wrap_mechss(sys, usfs)
%% function [eqn, oper] = mess_wrap_mechss(sys, usfs)
%
% Input
%   sys             mechss(M,C,K,B,F,G,D) a continuous-time first-order sparse state-space model
%                   object of the following form:
%                   M*x''(t) C*x'(t) + K*x(t) = B*u(t)
%                                        y(t) = F*x(t) + G*x'(t) + D*u(t)
%
%   usfs            string: name of folder containing the function handle set
%                   (optional, defaults to 'so_1')
%
% Output
%   eqn             struct contains data for equations
%                   M*x"(t) + E x'(t) + K*x(t)= B2*u(t)
%                                         y(t)= Cp*x(t) + Cv*x'(t) + D*u(t)
%
%                   eqn.M_ = M
%                   eqn.E_ = C
%                   eqn.K_ = K
%                   eqn.C  = |Cp Cv|
%                   eqn.D  = D
%
%                   for usfs = 'so_1':
%                           | 0  |
%                   eqn.B = | B2 |
%
%                   for usfs = 'so_2':
%                           | B2 |
%                   eqn.B = | 0  |
%
%
%   oper            struct contains function handles for operation with A and E
%
%

%% set oper
if not(exist('usfs', 'var'))
    oper = operatormanager('so_1');  %default setting
else
    oper = operatormanager(usfs);
end

%% set eqn
eqn.sys = sys;

if not(isempty(find(eqn.sys.M, 1)))
    eqn.M_ = eqn.sys.M;
end
if not(isempty(find(eqn.sys.C, 1)))
    eqn.E_ = eqn.sys.C;
end
if not(isempty(find(eqn.sys.K, 1)))
    eqn.K_ = eqn.sys.K;
end
if not(isempty(find(eqn.sys.B, 1)))
    switch(usfs)
        case 'so_1'
            eqn.B = [zeros(size(eqn.sys.B)); full(eqn.sys.B)];
        case  'so_2'
            eqn.B = [full(eqn.sys.B); zeros(size(eqn.sys.B))];
        otherwise
            warning('MESS:warning_arguments','eqn.B is only set for function handles ''so1'' and ''so2''');
    end
end

if not(isempty(eqn.sys.F))&& not(isempty(eqn.sys.G))
    eqn.C = [full(eqn.sys.F), full(eqn.sys.G)];
elseif not(isempty(eqn.sys.F))&& isempty(eqn.sys.G)
    eqn.C = [full(eqn.sys.F), zeros(size(eqn.sys.F))];
elseif isempty(eqn.sys.F)&& not(isempty(eqn.sys.G))
    eqn.C = [zeros(size(eqn.sys.G)), full(eqn.sys.G),];
else
    warning('MESS:warning_arguments','Neither Cp nor Cv is given. eqn.C will be missing');
end

if not(isempty(find(eqn.sys.D, 1)))
    eqn.D = eqn.sys.D;
end
