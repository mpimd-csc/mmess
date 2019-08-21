function [out, eqn, opts, oper] = adaptive_SDIRK43(eqn, opts, oper, h, L, t0)
% Solve the system E(t)' \dot{z}(t) = A(t)' z(t) with z(0) = L over 
% the interval [t0, t0 + h]. If t0 is omitted it is assumed to be 0.
% If eqn.type == 'N', instead solve E(t) \dot{z}(t) = A(t) z(t) .
% The method coefficients are taken from [1, p.100].
% 
% [1] Hairer, E. and Wanner, G., Solving Ordinary Differential Equations
% II: Stiff and Differential-Algebraic Problems, 2nd Ed., Springer, 2002.
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

[n, p] = size(L);
normL = norm(L);
errest_order = 3;
tend = t0 + h;
TOL = opts.exp_action.tol;
%     errest_order = errest_order - 1; % EPUS
%     TOL = tend * TOL;

c = [1/4, 3/4, 11/20, 1/2, 1];
aij= [     1/4,         0,      0,      0,   0;
    1/2,       1/4,      0,      0,   0;
    17/50,     -1/25,    1/4,      0,   0;
    371/1360, -137/2720, 15/544,    1/4,   0;
    25/24,    -49/48, 125/16, -85/12, 1/4];
beta = [25/24, -49/48, 125/16, -85/12, 1/4];
beta_err = [-3/16, -27/32, 25/32, 0, 1/4];

onep = ones(1, p);
B = zeros(5*p, 2*p);
for j = 1:5
    B((j-1)*p+1:j*p, 1:p) =  diag(beta(j)*onep);
    B((j-1)*p+1:j*p, p+1:2*p) =  diag(beta_err(j)*onep);
end

C = zeros(5*p, 5*p);
for j = 1:5
    for i = 1:j-1
        C((i-1)*p+1:i*p, (j-1)*p+1:j*p) = diag(aij(j,i)*onep);
    end
    C((j-1)*p+1:j*p, 1:p) =  diag(beta(j)*onep);
end

t = t0;
hj = (tend-t0) / 100; % Initial step size, guessing

if hj == 0
    out.Z = L;
    return
end
j = 1; % internal step
u = L;
while t(end) + hj <= tend
    K = zeros(n, 5*p);
    
    for i = 1:5
        [eqn, opts, oper] = ...
            opts.splitting.eval_matrix_functions(eqn, opts, oper, ...
                                                 t(end)+c(i)*hj);
        RHS = oper.mul_A(eqn, opts, eqn.type, ...
            (u + hj*K(:, 1:(i-1)*p)*C(1:(i-1)*p, (i-1)*p+1:i*p)), 'N');
        % Want to do (M - aij*h*A)\RHS, but only have support for
        % (A + pE)\RHS, so scale by -1/(aih*h)
        coeff = -1/(aij(i,i)*hj);
        RHS = coeff * RHS;
        
        K(:, (i-1)*p+1:i*p) = oper.sol_ApE(eqn, opts, eqn.type, ...
                                           coeff, eqn.type, RHS, 'N');
    end
    
    % 4th-order approximation:
    %     v2 = u + hj*(   25/24*k1 - 49/48*k2 +
    %                  + 125/16*k3 - 85/12*k4 + 1/4*k5);
    v = u + hj*K*B(:,1:p);
    % 3rd-order approximation:
    %     v = u + h*(   59/48*k1 - 17/96*k2 +
    %                + 225/32*k3 - 85/12*k4);
    % Difference between 4th- and 3rd-order approximations,
    % evaluated in a better way:
    %     errest = hj*norm(  -3/16*k1 -27/32*k2 +
    %                      + 25/32*k3 + 0*k4 + 1/4*k5);
    errest = hj/(1+normL)*norm(K*B(:,p+1:2*p));
    %         errest = errest / hj; % EPUS
    
    if errest > TOL % redo step
        hj = (0.9*TOL/errest)^(1/(errest_order)) * hj;
    else
        t(end+1) = t(end) + hj; %#ok<AGROW>
        
        % New step size
        hj = (0.95*TOL/errest)^(1/errest_order) * hj;
        
        j = j+1;
        u = v;
        
        if abs(t(end) - tend) < eps % We are at the final time
            out.Z = v;
            break
        end
        
        % Ensure that we end up precisely at tend
        if t(end) + hj > tend
            hj = tend-t(end);
        end
    end
end

end