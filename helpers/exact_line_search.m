function [ lambda ] = exact_line_search( W_old, DeltaK_old, W, DeltaK, S, S_old )
% Compute lambda for exact line search

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

%% Check input
if isempty(DeltaK_old)
    DeltaK_old = 0;
end
%% Compute scalar values
if isempty(S) && isempty(S_old)
    alpha = sum(sum(( (W_old' * W_old) ).^2, 1), 2) ...
        + sum(sum((DeltaK_old' * DeltaK_old).^2, 1), 2) ...
        - 2 * sum(sum((DeltaK_old' * W_old).^2, 1), 2);
    beta = sum(sum((W' * W).^2, 1), 2);
    delta = sum(sum((DeltaK' * DeltaK).^2, 1), 2);
    gamma = sum(sum((W_old' * W).^2, 1), 2) ...
        - sum(sum((DeltaK_old' * W).^2, 1), 2);
    epsilon = sum(sum((W_old' * DeltaK).^2, 1), 2) ...
        - sum(sum((DeltaK_old' * DeltaK).^2, 1), 2);
    zeta = sum(sum((DeltaK' * W).^2, 1), 2);
elseif not(isempty(S)) && not(isempty(S_old))
    S = diag(sqrt(S));
    S_old = diag(sqrt(S_old));
    alpha = sum(sum((S_old * (W_old' * W_old) * S_old).^2, 1), 2) ...
        + sum(sum((DeltaK_old' * DeltaK_old).^2, 1), 2) ...
        - 2 * sum(sum((DeltaK_old' * W_old * S_old).^2, 1), 2);
    beta = sum(sum((S * (W' * W) * S).^2, 1), 2);
    delta = sum(sum((DeltaK' * DeltaK).^2, 1), 2);
    gamma = sum(sum((S_old * W_old' * W * S).^2, 1), 2) ...
        - sum(sum((DeltaK_old' * W * S).^2, 1), 2);
    epsilon = sum(sum((S_old * W_old' * DeltaK).^2, 1), 2) ...
        - sum(sum((DeltaK_old' * DeltaK).^2, 1), 2);
    zeta = sum(sum((DeltaK' * W * S).^2, 1), 2);
else
    error('MESS:line_search', ...
        'Incorrect data for S and S_old.');
end

%% Compute lambda via eigenproblem
a0 = 2 * (gamma - alpha);
a1 = 2 * (alpha + beta - 2 * (gamma + epsilon));
a2 = 6 * (epsilon - zeta);
a3 = 4 * delta;
a = [a0 a1 a2 a3];
a = a / norm(a);

A = [0,     1,      0
    0,     0,      1
    -a(1),   -a(2),    -a(3)];
B = eye(3);
B(3, 3) = a(4);

% Octave does not support eig(A, B, 'qz') use eig(A,B) as fallback in this case.
try 
  lambda = eig(A, B, 'qz');
catch
  lambda = eig(A, B);
end

lambda = lambda(not(imag(lambda)));
lambda = lambda(lambda >= 0);
lambda = lambda(lambda <= 2);

if isempty(lambda)
    lambda = 0;
    warning('MESS:exact_line_search', ...
                    'Could not find a stepsize lambda.');
elseif size(lambda, 1) > 1
    f = @(t) ((1 - t).^2) * alpha + (t.^2) * beta + (t.^4) * delta ...
        + (2 * t .* (1 - t)) * gamma - (2 * t.^2 .* (1 - t)) * epsilon ...
        - (2 * t.^3) * zeta;
    [~, I] = min(f(lambda));
    lambda = lambda(I);
end
%% With fminbnd
% f = @(t) ((1 - t).^2) * alpha + (t.^2) * beta + (t.^4) * delta ...
%     + (2 * t .* (1 - t)) * gamma - (2 * t.^2 .* (1 - t)) * epsilon ...
%     - (2 * t.^3) * zeta;
% tol = 1e-12;
% lambda = fminbnd(@(t) f(t), 0, 2, optimset('TolX', tol));

%% Test as eigenproblem (DOI. 10.1137/S0895479899365720)

% A1 = [0,     1,      0
%     0,     0,      1
%     -a0,   -a1,    -a2];
% B1 = eye(3);
% B1(3, 3) = a3;
% E1 = eig(A1, B1); % inaccurate
% E1 = E1(not(imag(E1)));
% E1 = E1(E1 >= 0);
% E1 = E1(E1 <= 2);

% E2 = eig(A1, B1, 'qz'); % inaccurate
% E2 = E2(not(imag(E2)));
% E2 = E2(E2 >= 0);
% E2 = E2(E2 <= 2);
 
% C1 = [0,         1,          0
%     0,          0,          1
%     -a0 / a3,   -a1 / a3,   -a2 / a3];
% E3 = eig(C1);
% E3 = E3(not(imag(E3)));
% E3 = E3(E3 >= 0);
% E3 = E3(E3 <= 2);


% E4 = eig(A2, B2);
% E4 = E4(not(imag(E4)));
% E4 = E4(E4 >= 0);
% E4 = E4(E4 <= 2);

% 
% C2 = [0,            1,              0
%     0,              0,              1
%     -a(1) / a(4),   -a(2) / a(4),   -a(3) / a(4)];
% E6 = eig(C2); % seems to be most reliable in tests
% E6 = E6(not(imag(E6)));
% E6 = E6(E6 >= 0);
% E6 = E6(E6 <= 2);

end

