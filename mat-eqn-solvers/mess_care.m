function [Z, D, K] = mess_care(A, B, C, S, E)
% mess_care Solve continuous-time Riccati equations with 
%           sparse coefficients.
% 
%   [Z, K] = mess_care(A, B, C) solves the Riccati matrix equation:
% 
%        A'*Z*Z' + Z*Z'*A - Z*Z'*B*B'*Z*Z' + C'*C = 0
%
%        K is the feedback matrix K = B'*Z*Z'
%        To omit the computation of Z use: 
%   K = mess_care(A, B, C)
%        To get only the solution factor Z as output use:
%   [Z, ~] = mess_care(A, B, C)
% 
%
%   [Z, K] = mess_care(A, B, C, [], E) solves the generalized Riccati 
%        equation:
% 
%        A'*Z*Z'*E + E'*Z*Z'*A - E'*Z*Z'*B*B'*Z*Z'*E + C'*C = 0 
%
%        K is the feedback matrix K = B'*Z*Z'*E
%        To omit the computation of Z use: 
%   K = mess_care(A, B, C, [], E)
%        To get only the solution factor Z as output use:
%   [Z, ~] = mess_care(A, B, C, [], E)
% 
%
%   [Z, D, K] = mess_care(A, B, C, S) solves the Riccati matrix equation
%       in ZDZ^T formulation:
% 
%        A'*Z*D*Z' + Z*D*Z'*A - Z*D*Z'*B*B'*Z*D*Z' + C'*S*C = 0
%
%        K is the feedback matrix K = B'*Z*D*Z'
%        To omit the computation of Z and D use: 
%   K = mess_care(A, B, C, S)
%        To get only the solution factors Z and D as output use:
%   [Z, D] = mess_care(A, B, C, S)
% 
%
%   [Z, D, K] = mess_care(A, B, C, S, E) solves the generalized Riccati 
%       equation in ZDZ^T formulation:
% 
%        A'*Z*D*Z'*E + E'*Z*D*Z'*A - E'*Z*D*Z'*B*B'*Z*D*Z'*E + C'*S*C = 0 
%
%        K is the feedback matrix K = B'*Z*D*Z'*E
%        To omit the computation of Z and D use: 
%   K = mess_care(A, B, C, S, E)
%        To get only the solution factor Z as output use:
%   [Z, D] = mess_care(A, B, C, S, E)
% 
%
%   For the dense case see also care, dare.

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
oper = operatormanager('default');

%% Options
ni = nargin;
no = nargout;
if ni < 4
    S = [];
end
% ADI
opts.adi.info = 0;
opts.adi.res_tol = 1e-12;
opts.adi.rel_diff_tol = 0;
opts.adi.maxiter = 100;
n = size(A, 1);
opts.shifts.num_Ritz = min(n - 2, 50);
opts.shifts.num_hRitz = min(n - 2, 25);
if no == 1
    opts.shifts.method = 'heur';
    opts.shifts.num_desired = min(floor((opts.shifts.num_Ritz + opts.shifts.num_hRitz) / 2) - 2, 25);
else
    opts.shifts.method = 'projection';
    opts.shifts.num_desired = max(6, size(B, 2));
end
opts.norm = 'fro';
opts.adi.compute_sol_fac = 1;
opts.adi.accumulateK = 1;
opts.adi.accumulateDeltaK = 0;
if no == 1
    opts.adi.compute_sol_fac = 0;
else
    opts.adi.compute_sol_fac = 1;
end
% NM
opts.nm.maxiter = 25;
opts.nm.res_tol = 1e-11;
opts.nm.rel_diff_tol = 0;
opts.nm.info = 0;
opts.nm.accumulateRes = 1;
opts.nm.linesearch = 1;
opts.norm = 'fro';


%% Equation type
eqn.type = 'T';
if isempty(S)
    eqn.A_ = A;
    eqn.B = B;
    eqn.C = C;
    if ni == 3
        eqn.haveE = 0;
    elseif ni == 5
        eqn.haveE = 1;
        eqn.E_ = E;
    else
        error('MESS:notimplemented', 'Wrong number of input arguments');
    end
else % ZDZ^T case 
    opts.LDL_T = 1;
    eqn.A_ = A;
    eqn.B = B;
    eqn.C = C;
    eqn.S = S;
    if ni == 4
        eqn.haveE = 0;
    elseif ni == 5
        eqn.haveE = 1;
        eqn.E_ = E;
    else
        error('MESS:notimplemented', 'Feature not yet implemented!');
    end
end

%% Solve Equation
out = mess_lrnm(eqn, opts, oper);

%% Prepare output
if no >= 2 
    Z = out.Z;
end
if (not(isempty(S))) && (no >= 2) % ZDZ^T case 
    D = out.D;
    if no == 3
       K = out.K; 
    end
elseif no == 2
    D = out.K;
elseif no == 1
    Z = out.K;
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Comments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Projection shifts require Z to be computed. They use the last columns in
% Z to project down the matrices.
% With only one output argument, Z should not be computed 
% (opts.adi.compute_sol_fac = 0). Thus, we use heurisic shifts. 
