function [Z, D] = mess_lyap(varargin)
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
%    If S is emptry, matrices A,B and E can be given as Z = mess_lyap(sys)
%    with sys = sparss(A, B , C_ ,D , E) a continuous-time first-order sparse
%    state-space object of the following form:
%                   E*x'(t) = A*x(t)  + B*u(t)
%                   y(t)    = C_*x(t) + D*u(t)
%
%    For the dense case see also lyap, lyapchol, dlyap.

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2021 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%


%% Usfs
oper=operatormanager('default');


%% Options
ni = nargin;
no = nargout;

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

%% Decide if sparss or single matrices were passed in
if ni == 1 && isa(varargin{1}, 'sparss')
    [eqn, oper] = mess_wrap_sparss(varargin{1});
    % set eqn properties
    eqn.type = 'N';
    if(exist('eqn.C', 'var'))
        warning('MESS:ignored',...
                'C is supposed to be empty. Data is ignored.');
    end
    if(exist('eqn.D', 'var'))
        warning('MESS:ignored',...
                'D is supposed to be empty. Data is ignored.');
    end
    eqn.G = eqn.B;
else
    if ni < 5
        S = [];
    else
        S = varargin{4};
    end
    %% Equation type
    eqn.type = 'N';
    if no == 1
        if not(isempty(S))
            warning('MESS:ignored',...
                'Fourth argument is supposed to be empty. Data is ignored.');
        end
        eqn.A_ = varargin{1};
        eqn.G = varargin{2};
        if ni == 2
            eqn.haveE = 0;
        elseif ni == 5
            if  not(isempty(varargin{3}))
                warning('MESS:ignored',...
                    'Third argument is supposed to be empty. Data is ignored.');
            end
            eqn.haveE = 1;
            eqn.E_ = varargin{5};
        else
            error('MESS:notimplemented', 'Feature not yet implemented!');
        end
    elseif no == 2 % ZDZ^T case
        opts.LDL_T = 1;
        eqn.A_ = varargin{1};
        eqn.G = varargin{2};
        eqn.S = varargin{4};
        if ni == 4
            if not(isempty(varargin{3}))
                warning('MESS:ignored',...
                    'Third argument is supposed to be empty. Data is ignored.');
            end
            eqn.haveE = 0;
        elseif ni == 5
            if not(isempty(varargin{3}))
                warning('MESS:ignored',...
                    'Third argument is supposed to be empty. Data is ignored.');
            end
            eqn.haveE = 1;
            eqn.E_ = varargin{5};
        else
            error('MESS:notimplemented', 'Feature not yet implemented!');
        end
    else
        error('MESS:notimplemented', 'Feature not yet implemented!');
    end
    eqn.B = varargin{2};
end


%% Solve Equation
out = mess_lradi(eqn, opts, oper);
Z   = out.Z;

%% Prepare output
if no == 2 % ZDZ^T case
    D = out.D;
end
