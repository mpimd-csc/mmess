function [Z, D] = mess_lyap(A, B, C, S, E)
% mess_lyap Solve continuous-time Lyapunov equations with 
%           sparse coefficients.
% 
%    Z = mess_lyap(A, B) solves the Lyapunov matrix equation:
% 
%        A*Z*Z' + Z*Z'*A' + B*B' = 0
% 
%    [Z, D, Y] = mess_lyap(A, B, C) solves the Sylvester equation:
% 
%        A*Z*D*Z' + Z*D*Y'*B + C = 0 (NOT YET IMPLEMENTED)
% 
%    Z = mess_lyap(A, B, [], [], E) solves the generalized Lyapunov 
%        equation:
% 
%        A*Z*Z'*E' + E*Z*Z'*A' + B*B' = 0 
% 
%    [Z, D] = mess_lyap(A, B, [], S) solves the Lyapunov matrix equation
%       in ZDZ^T formulation:
% 
%        A*Z*D*Z' + Z*D*Z'*A' + B*S*B' = 0
% 
%    [Z, D] = mess_lyap(A, B, [], S, E) solves the generalized Lyapunov 
%       equation in ZDZ^T formulation:
% 
%        A*Z*D*Z'*E' + E*Z*D*Z'*A' + B*S*B' = 0 
% 
%    For the dense case see also lyap, lyapchol, dlyap.

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

%% Usfs
oper=operatormanager('default');

%% Options
ni=nargin;
no = nargout;
if ni < 4
    S = [];
end
opts.adi.info=0;
opts.adi.res_tol=1e-12;
opts.adi.rel_diff_tol=0;
opts.adi.maxiter=100;
opts.shifts.num_Ritz=50;
opts.shifts.num_hRitz=25;
opts.shifts.method = 'projection';
opts.shifts.num_desired = 6;
opts.norm = 'fro';
opts.adi.compute_sol_fac = 1;


%% Equation type
eqn.type = 'N';
if no == 1
    if not(isempty(S))
        warning('MESS:ignored',...
            'Fourth argument is supposed to be empty. Data is ignored.');
    end
    eqn.A_ = A;
    eqn.G = B;
    if ni == 2
        eqn.haveE = 0;
    elseif ni == 5
        if not(isempty(C))
            warning('MESS:ignored',...
                'Third argument is supposed to be empty. Data is ignored.');
        end
        eqn.haveE = 1;
        eqn.E_ = E;
    else
        error('MESS:notimplemented', 'Feature not yet implemented!');
    end
elseif no == 2 % ZDZ^T case
    opts.LDL_T = 1;
    eqn.A_ = A;
    eqn.G = B;
    eqn.S = S;
    if ni == 4
        if not(isempty(C))
            warning('MESS:ignored',...
                'Third argument is supposed to be empty. Data is ignored.');
        end
        eqn.haveE = 0;
    elseif ni == 5
        if not(isempty(C))
            warning('MESS:ignored',...
                'Third argument is supposed to be empty. Data is ignored.');
        end
        eqn.haveE = 1;
        eqn.E_ = E;
    else
        error('MESS:notimplemented', 'Feature not yet implemented!');
    end
else
    error('MESS:notimplemented', 'Feature not yet implemented!');
end
eqn.B = B;

%% Shift parameter
opts.shifts.p = mess_para(eqn, opts, oper);

%% Solve Equation
out = mess_lradi(eqn, opts, oper);
Z   = out.Z;

%% Prepare output
if no == 2 % ZDZ^T case 
    D = out.D;
end
