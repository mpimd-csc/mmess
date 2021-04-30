function [eqn, oper] = mess_wrap_sparss(sys, usfs)
%% function [eqn, oper] = mess_wrap_sparss(sys, usfs)
%
% Input
%   sys             sys = sparss(A,B,C,D,E) a continuous-time first-order
%                   sparse state-space model object of the following form:
%                   E*x'(t) = A*x(t) + B*u(t)
%                   y(t)    = C*x(t) + D*u(t)
%
%   usfs            string: name of folder containing the function handle set
%                   (optional, defaults to 'default')
%
% Output
%   eqn             struct contains data for equations
%
%   oper            struct contains function handles for operation with A and E
%
%


%% set oper
if not(exist('usfs', 'var'))
    oper = operatormanager('default');  %default setting
else
    oper = operatormanager(usfs);
end

%% set eqn
eqn.sys = sys;

if not(isempty(find(eqn.sys.A, 1)))
    eqn.A_ = eqn.sys.A;
end
if not(isempty(find(eqn.sys.B, 1)))
    eqn.B = full(eqn.sys.B);
end
if not(isempty(find(eqn.sys.C, 1)))
    eqn.C = full(eqn.sys.C);
end
if not(isempty(find(eqn.sys.D, 1)))
    eqn.D = eqn.sys.D;
end
if not(isempty(find(eqn.sys.E, 1)))
    eqn.E_ = eqn.sys.E;
    eqn.haveE = 1;
else
    eqn.haveE = 0;
end

