function [ROM, outinfo, eqn, opts, oper] = ...
    mess_balanced_truncation_bilinear(eqn, opts, oper)
% Lyapunov Balanced truncation for descriptor systems with invertible E.
%
%  [out, eqn, opts, oper] = mess_balanced_truncation(eqn, opts, oper)
%
% Input
%   eqn                 struct contains data for equations
%
%   opts                struct contains parameters for the algorithm
%
%   oper                struct contains function handles for operation
%                       with A, E and N
%
% Output
%   out                 struct containing output information
%
%   eqn                 struct contains data for equations
%
%   opts                struct contains parameters for the algorithm
%
%   oper                struct contains function handles for operation
%                       with A, E and N
%
%
% Input fields in struct eqn:
%   eqn.B       dense (n x m1) matrix B
%
%   eqn.C       dense (m2 x n) matrix C
%
%   Depending on the operator chosen by the operatormanager, additional
%   fields may be needed. For the "default", e.g., eqn.A_ and eqn.E_ hold
%   the A and E matrices. Moreover for the bilinear models we need
%   eqn.N_      Cell with N_k = N{k} for k = 1,2, ...
%               dense (n x n)  matrix N_k
%               (if all N_k are given in one large matrix it will be
%                transformed to a cell. output is matrix again)
%
%
% Input fields in struct opts:
%
%  opts       options structure that can be used to pass setting to the
%             LRADI, ADI shift computation, or the square root method (optional)
%             (see corresponding routines for additional information)
%
%
% Output fields in struct out:
%   Er, Ar, Br, Cr, Nr_          the reduced order model matrices
%                                (out.Nr_ is given as cell array)
%
%
%   out.outB_lyapunov_bilinear   output information of the lyapunov solver
%                                A*Z*Z'*E' + E*Z*Z'*A' + Sum_N_k*Z*N_k' + B*B' = 0 (N - Case)
%                                (see corresponding routine for additional information)
%
%
%   out.outC_lyapunov_bilinear   output information of the lyapunov solver
%                                A'*Z*Z'*E + E'*Z*Z'*A + Sum_N_k'*Z*N_k + C'*C = 0 (T - Case)
%                                (see corresponding routine for additional information)
%
%
%   out.TL and out.TR            left and right truncation matrices
%   out.hsv                      computed Hankel singular values
%                                (by the square_root_method)
%
% References:
%
% [1] P. Benner, P. Goyal, Balanced Truncation Model Order Reduction For
%     Quadratic-Bilinear Control Systems, e-prints 1705.00160, arXiv, math.OC
%     (2017). URL https://arxiv.org/abs/1705.00160


%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%


%% Check oper and initialize parameters
% operations are done by the default set of user supplied functions
if nargin < 3
    oper = operatormanager('default');
end
% make sure we use default usfs as none of the others supports mul_N so far
if not(isequal(oper.name, 'default'))
    error('MESS:notimplemented', ...
        [oper.name, ' usfs are not supported in this function.']);
end

% Initialize variables
n = oper.size(eqn, opts);

% BT tolerance and maximum order for the ROM
if not(isfield(opts,'srm')) || not(isfield(opts.srm, 'tol'))
    opts.srm.tol = 1e-5;
end

if not(isfield(opts.srm, 'max_ord'))
    opts.srm.max_ord = n;
end

if not(isfield(opts.srm, 'info'))
    opts.srm.info = 0;
end

% some control settings for the LRADI
if not(isfield(opts,'adi')) || not(isfield(opts.adi, 'maxiter'))
    opts.adi.maxiter = 100;
end

if not(isfield(opts.adi, 'res_tol'))
    opts.adi.res_tol = 1e-9;
end

if not(isfield(opts.adi, 'rel_diff_tol'))
    opts.adi.rel_diff_tol = 1e-16;
end

if not(isfield(opts, 'norm'))
    opts.norm = 'fro';
end

%% Check Problem data
if not(isfield(eqn,'haveE'))
    eqn.haveE = 0;
    warning('MESS:control_data', ...
            ['Missing or corrupted eqn.haveE field.', ...
             'Switching to default: 0']);
end

% setup USFS
[eqn, opts, oper] = oper.mul_A_pre(eqn, opts, oper);
[eqn, opts, oper] = oper.mul_N_pre(eqn, opts, oper);


%% make sure proper shift selection parameters are given
% If not set outside, we use projection shifts
if not(isfield(opts,'shifts')) || not(isfield(opts.shifts, 'method'))
    opts.shifts.method = 'projection';
end

if not(isfield(opts.shifts, 'num_desired'))
    opts.shifts.num_desired = max(5, min(size(eqn.B, 2), size(eqn.C, 1)));
end

if not(isfield(opts.shifts, 'b0'))
    opts.shifts.b0 = ones(n, 1);
end


%% Truncated controllability Gramian
eqn.type = 'N';

[outB_lyapunov_bilinear, eqn, opts, oper] = ...
    mess_lyapunov_bilinear(eqn, opts, oper);

outinfo.outB_lyapunov_bilinear = outB_lyapunov_bilinear;

%% Truncated observability Gramian
eqn.type = 'T';

[outC_lyapunov_bilinear, eqn, opts, oper] = ...
    mess_lyapunov_bilinear(eqn, opts, oper);

outinfo.outC_lyapunov_bilinear = outC_lyapunov_bilinear;

%% Square root method
[outinfo.TL ,outinfo.TR, outinfo.hsv] = mess_square_root_method(eqn, opts , ...
    oper, outB_lyapunov_bilinear.Z, outC_lyapunov_bilinear.Z);

%% compute ROM matrices
ROM.A = outinfo.TL' * oper.mul_A(eqn, opts, 'N', outinfo.TR, 'N');
ROM.B = outinfo.TL' * eqn.B;
ROM.C = eqn.C * outinfo.TR;
ROM.E = eye(size(ROM.A, 1));

numberOf_N_matrices = length(eqn.N_);

for currentN_k = 1 : numberOf_N_matrices

    ROM.N{currentN_k} = outinfo.TL' * oper.mul_N(eqn, opts, 'N', ...
        outinfo.TR, 'N', currentN_k);
end

%% clean up usfs
[eqn, opts, oper] = oper.mul_A_post(eqn, opts, oper);
[eqn, opts, oper] = oper.mul_N_post(eqn, opts, oper);

