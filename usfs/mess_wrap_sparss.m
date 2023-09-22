function [eqn, opts, oper] = mess_wrap_sparss(sys, opts, usfs)
%% function [eqn, opts, oper] = mess_wrap_sparss(sys, opts, usfs)
%
% Input
%   sys             sys = sparss(A,B,C,D,E) a continuous-time first-order
%                   sparse state-space model object of the following form:
%                   E*x'(t) = A*x(t) + B*u(t)
%                   y(t)    = C*x(t) + D*u(t)
%
%   opts            transit argument required by logger functions
%
%   usfs            string: name of folder containing the function handle set
%                   (optional, defaults to 'default')
%
% Output
%   eqn             struct contains data for equations
%
%   opts            transit argument required by logger functions
%
%   oper            struct contains function handles for operation with A and E
%
%

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%%
narginchk(2, 3);

%% set oper
if nargin < 3 || isempty(usfs)
    [oper, opts] = operatormanager(opts, 'default');  % default setting
else
    [oper, opts] = operatormanager(opts, usfs);
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
    eqn.haveE = true;
else
    eqn.haveE = false;
end
