function [Z, D, K] = mess_care(varargin)
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
%   If S is empty, matrices A,B and E can be given as Z = mess_lyap(sys)
%   with sys = sparss(A, B , C_ ,D , E) a continuous-time first-order sparse
%   state-space object of the following form:
%                   E*x'(t) = A*x(t)  + B*u(t)
%                   y(t)    = C_*x(t) + D*u(t)
%
%   For the dense case see also care, dare.

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%


%% Usfs
oper = operatormanager('default');

%% Options
ni = nargin;
no = nargout;

%% Equation type
if (ni == 1) && isa(varargin{1}, 'sparss')
    [eqn, oper] = mess_wrap_sparss(varargin{1});
    eqn.type = 'T';
    if(exist('eqn.D', 'var'))
        warning('MESS:ignored',...
                'D is supposed to be empty. Data is ignored.');
    end
    S = [];
else
    eqn.type = 'T';
    if ni < 4
        S = [];
    else
        S = varargin{4};
    end
    if isempty(S) % Z*Z' case.
        eqn.A_ = varargin{1};
        eqn.B  = varargin{2};
        eqn.C  = varargin{3};

        if ni == 3
            eqn.haveE = 0;
        elseif ni == 5
            eqn.haveE = 1;
            eqn.E_ = varargin{5};
        else
            error('MESS:notimplemented', 'Wrong number of input arguments');
        end
    else % Z*D*Z' case.
        opts.LDL_T = 1;
        eqn.A_ = varargin{1};
        eqn.B = varargin{2};
        eqn.C = varargin{3};
        eqn.S = varargin{4};

        if ni == 4
            eqn.haveE = 0;
        elseif ni == 5
            eqn.haveE = 1;
            eqn.E_ = varargin{5};
        else
            error('MESS:notimplemented', 'Feature not yet implemented!');
        end
    end
end

% Global.
opts.norm = 'fro';

% Shifts.
n  = size(eqn.A_, 1);
s1 = min(n - 2, 50);
s2 = min(n - 2, 25);

opts.shifts.num_desired       = min(floor((s1 + s2) / 2) - 2, 25);
opts.shifts.history           = opts.shifts.num_desired * size(eqn.C, 1);
opts.shifts.method            = 'gen-ham-opti';
opts.shifts.naive_update_mode = false;

% RADI.
opts.radi.maxiter      = 300;
opts.radi.res_tol      = 1.0e-11;
opts.radi.rel_diff_tol = 0;
opts.radi.info         = 0;
opts.radi.trunc_tol    = 1.0e-13;

switch no

    case 1
        % Compute only K.
        opts.radi.compute_sol_fac = 0;
        opts.radi.get_ZZt         = 0;
    case 2
        % Compute K and Z in Z*Z' format.
        opts.radi.compute_sol_fac = 1;
        opts.radi.get_ZZt         = 1;
    otherwise
        % Compute K, Z and D in Z*D*Z' format.
        opts.radi.compute_sol_fac = 1;
        opts.radi.get_ZZt         = 0;
end

%% Solve Equation
out = mess_lrradi(eqn, opts, oper);

if out.res(end) > opts.radi.res_tol
    warning('MESS:convergence', ...
            ['Convergence of solution only up to relative residual of %e!\n' ...
             'Check mess_lrnm and mess_lrradi for customizable solvers.'], ...
            out.res(end));
end

%% Prepare output
if no >= 2
    Z = out.Z;
end

if (not(isempty(S))) && (no >= 2) % Z*D*Z' case.
    D = out.D;
    if no == 3
       K = out.K;
    end
elseif no == 2 % Z*Z' case and K.
    D = out.K;
elseif no == 1 % Only K.
    Z = out.K;
end
