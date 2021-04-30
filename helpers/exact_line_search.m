function [ lambda ] = exact_line_search( W_old, DeltaK_old, W, DeltaK, S, S_old )
% Compute lambda for exact line search

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2021 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
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

end
%Here we keep alternative formulations that we tested in case this routine
%is investigated again in the future.
%% This is the original code by Heiko Wichelt with fminbnd
% f = @(t) ((1 - t).^2) * alpha + (t.^2) * beta + (t.^4) * delta ...
%     + (2 * t .* (1 - t)) * gamma - (2 * t.^2 .* (1 - t)) * epsilon ...
%     - (2 * t.^3) * zeta;
% tol = 1e-12;
% lambda = fminbnd(@(t) f(t), 0, 2, optimset('TolX', tol));

%% As eigenproblem
% We tested different variants from (DOI. 10.1137/S0895479899365720) to
% formulate the step-size computation as an eigenvalue problem for accuracy
% and performance. Below is code from testing the variants. The tests were
% done in 2015 (MESS version < 1.0). We chose variant 5 for the
% implementation above since it produced a step-size closest to the one from
% the original approach in most of the test cases.

% %% Test as eigenproblem
% fprintf('\ntesting line search:\n');
% fprintf('alpha: %e\n', alpha);
% fprintf('beta: %e\n', beta);
% fprintf('\t\tdelta: %e\n', delta);
% fprintf('gamma: %e\n', gamma);
% fprintf('epsilon: %e\n', epsilon);
% fprintf('zeta: %e\n', zeta);
%
% a0 = 2 * (gamma - alpha);
% a1 = 2 * (alpha + beta - 2 * (gamma + epsilon));
% a2 = 6 * (epsilon - zeta);
% a3 = 4 * delta;
% fprintf('\na0: %e\n', a0);
% fprintf('a1: %e\n', a1);
% fprintf('a2: %e\n', a2);
% fprintf('a3: %e\n', a3);
%
% A1 = [0,     1,      0
%     0,     0,      1
%     -a0,   -a1,    -a2];
% B1 = eye(3);
% B1(3, 3) = a3;
% E1 = eig(A1, B1);
% E1 = E1(not(imag(E1)));
% E1 = E1(E1 >= 0);
% E1 = E1(E1 <= 2);
% for i = 1 : length(E1)
%     fprintf('1: diff: %e \t val diff: %e\n', abs(E1(i) - lambda), ...
%         f(lambda) - f(E1(i)));
% end
% E2 = eig(A1, B1, 'qz');
% E2 = E2(not(imag(E2)));
% E2 = E2(E2 >= 0);
% E2 = E2(E2 <= 2);
% for i = 1 : length(E2)
%     fprintf('2: diff: %e \t val diff: %e\n', abs(E2(i) - lambda), ...
%         f(lambda) - f(E2(i)));
% end
%
% C1 = [0,         1,          0
%     0,          0,          1
%     -a0 / a3,   -a1 / a3,   -a2 / a3];
% E3 = eig(C1);
% E3 = E3(not(imag(E3)));
% E3 = E3(E3 >= 0);
% E3 = E3(E3 <= 2);
% for i = 1 : length(E3)
%     fprintf('3: diff: %e \t val diff: %e\n', abs(E3(i) - lambda), ...
%         f(lambda) - f(E3(i)));
% end
%
% a = [a0 a1 a2 a3];
% a = a / norm(a);
%
% A2 = [0,     1,      0
%     0,     0,      1
%     -a(1),   -a(2),    -a(3)];
% B2 = eye(3);
% B2(3, 3) = a(4);
% E4 = eig(A2, B2);
% E4 = E4(not(imag(E4)));
% E4 = E4(E4 >= 0);
% E4 = E4(E4 <= 2);
% for i = 1 : length(E4)
%     fprintf('4: diff: %e \t val diff: %e\n', abs(E4(i) - lambda), ...
%         f(lambda) - f(E4(i)));
% end
% E5 = eig(A2, B2, 'qz');
% E5 = E5(not(imag(E5)));
% E5 = E5(E5 >= 0);
% E5 = E5(E5 <= 2);
% for i = 1 : length(E5)
%     fprintf('5: diff: %e \t val diff: %e\n', abs(E5(i) - lambda), ...
%         f(lambda) - f(E5(i)));
% end
%
% C2 = [0,            1,              0
%     0,              0,              1
%     -a(1) / a(4),   -a(2) / a(4),   -a(3) / a(4)];
% E6 = eig(C2);
% E6 = E6(not(imag(E6)));
% E6 = E6(E6 >= 0);
% E6 = E6(E6 <= 2);
% for i = 1 : length(E6)
%     fprintf('6: diff: %e \t val diff: %e\n', abs(E6(i) - lambda), ...
%         f(lambda) - f(E6(i)));
% end
