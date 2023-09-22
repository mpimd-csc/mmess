function y = lyapunov(Z, x, eqn, oper, opts, D)
% Computes matrix vector product with the Lyapunov operator.
%
% Input:
%  Z         Low-rank solution factor of the Riccati equation
%  x         vector for matrix vector product
%  eqn       structure with data for A, E and fields W
%  oper      structure contains function handles for operations with A, E
%  opts      full options structure (passed on to function handles in oper)
%
%  D         solution factor D for the LDL^T
%            decomposition, i.e., opts.LDL_T = true
%  eqn.type 'N' or 'T' for type of Lyapunov equation
%
% Output:
%  y        result of matrix vector product

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

if eqn.type == 'N'
    adjoint = 'T';
else
    adjoint = 'N';
end

if isempty(D) && opts.LDL_T
    mess_err(opts, 'lyapunov:LDL^T formulation needs D to get passed.');
end

if eqn.haveE
    if opts.LDL_T
        z = Z * mess_LDL_mul_D(eqn, D, Z' * oper.mul_E(eqn, opts, adjoint, x, 'N'));
    else
        z = Z * (Z' * (oper.mul_E(eqn, opts, adjoint, x, 'N')));
    end
else
    if opts.LDL_T
        z = Z * mess_LDL_mul_D(eqn, D, Z' * x);
    else
        z = Z * (Z' * x);
    end
end

y1 = oper.mul_A(eqn, opts, eqn.type, z, 'N');
y2 = oper.mul_A(eqn, opts, adjoint, x, 'N');

if eqn.haveUV
    if eqn.type == 'N'
        y1 = y1 + eqn.pU * (eqn.V' * z);
        y2 = y2 + eqn.V * (eqn.U' * x);
    else
        y1 = y1 + eqn.pV * (eqn.U' * z);
        y2 = y2 + eqn.U * (eqn.V' * x);
    end
end

if opts.LDL_T
    y2 = Z * mess_LDL_mul_D(eqn, D, Z' * y2);
else
    y2 = Z * (Z' * y2);
end

if eqn.haveE
    y2 = oper.mul_E(eqn, opts, eqn.type, y2, 'N');
end

y = y1 + y2;

if opts.LDL_T
    y = y + eqn.W * (eqn.T * (eqn.W' * x));
else
    y = y + eqn.W * (eqn.W' * x);
end

% in case of Rosenbrock we get a -1/(2*timestep)*(E'*Z*Z'*E)
% from both F and F'
if isfield(opts, 'rosenbrock') && not(isempty(opts.rosenbrock))
    if eqn.haveE       % generalized equations
        y = y - (1.0 / opts.rosenbrock.stepsize) * oper.mul_E(eqn, opts, eqn.type, z, 'N');
    else
        y = y - (1.0 / opts.rosenbrock.stepsize) * z;
    end
end
