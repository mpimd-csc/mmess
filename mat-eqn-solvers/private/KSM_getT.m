function [eqn, opts, oper] = KSM_getT(eqn, opts, oper)
% function [eqn,opts,oper] = KSM_getT(eqn,opts,oper)
% Function that computes the update the projection of A onto the current
% subspace
%
% Input and output:
%
%   eqn       structure containing equation data
%
%   opts      structure containing parameters for the algorithm
%
%   oper      contains function handles with operations for A and E
%
%   opts.KSM.compute_struct  structure that contain all the useful
%                                  information computed in the previous
%                                  iterations
%
%
% Output:
%
%  opts.KSM.compute_struct       structure containing the following:
%
%  opts.KSM.compute_struct.V     basis of the current space
%
%  opts.KSM.compute_struct.H     matrix collecting the coefficients
%                                      stemming from the Arnoldi
%                                      (or Lanczos) process
%
%  opts.KSM.compute_struct.T     projection of A onto the current
%                                      space, namely T=V'*A*V
%
%  opts.KSM.compute_struct.L,    matrices needed to recover the
%                                      columns of T
%  opts.KSM.compute_struct.rho   from the columns of H
%
%  opts.KSM.compute_struct.beta  p-by-p matrix such that
%                                      B=V_1*compute_struct.beta
%

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

[eqn, opts, oper] = oper.mul_A_pre(eqn, opts, oper);
p = size(eqn.W, 2);

VV = opts.KSM.compute_struct.V;
it = opts.KSM.compute_struct.it;

if strcmp('EK', opts.KSM.space)
    Lnew = zeros(2 * p * (it + 1), 2 * p);
    Tnew = Lnew;
elseif strcmp('RK', opts.KSM.space)
    Tnew = zeros(p * (it + 1), p);
end

if strcmp('EK', opts.KSM.space)

    if opts.KSM.explicit_proj
        T = opts.KSM.compute_struct.T;
        AVnew = oper.mul_A(eqn, opts, eqn.type, ...
                           VV(:, 2 * p * (it - 1) + 1:2 * p * it), 'N');
        Tnew = VV' * AVnew;

        % expand T
        if it == 1
            opts.KSM.compute_struct.T = Tnew;
        else
            opts.KSM.compute_struct.T = ...
                [[T; zeros(2 * p, 2 * p * (it - 1))], Tnew];
        end

    else

        L = opts.KSM.compute_struct.L;
        T = opts.KSM.compute_struct.T;
        H = opts.KSM.compute_struct.H;
        hinv = inv(H(2 * p * it + 1:2 * p * (it + 1), :));

        % Compute the projected matrix T=V'*A*V from the columns of H
        if it == 1
            ibeta = inv(opts.KSM.compute_struct.beta);
            Lnew = [H(1:3 * p, 1:p) / ibeta(1:p, 1:p), ...
                    speye(3 * p, p) / ibeta(1:p, 1:p)] * ...
                   ibeta(1:2 * p, p + 1:2 * p);
        else
            Lnew(1:2 * p * it, 1:p) = L(1:2 * p * it, p + 1:2 * p) + ...
                H(1:2 * p * it, 1:p) * opts.KSM.compute_struct.rho;
            Lnew(2 * p * it + 1:2 * p * (it + 1), 1:p) = ...
                H(2 * p * it + 1:2 * p * (it + 1), 1:p) * ...
                opts.KSM.compute_struct.rho;
        end

        % the odd columns of T correspond to the odd column of H
        Tnew(1:2 * p * (it + 1), 1:p) = H(1:2 * p * (it + 1), 1:p); % odd columns

        % the even columns of T correspond to the even column of L
        % notice that the last pX2p block of such columns is a zero
        Tnew(1:2 * p * it + p, p + 1:2 * p) = Lnew(1:2 * p * it + p, 1:p);

        % expand T
        if it == 1
            opts.KSM.compute_struct.T = Tnew;
        else
            opts.KSM.compute_struct.T = ...
                [[T; zeros(2 * p, 2 * p * (it - 1))], Tnew];
        end

        I = eye(2 * p * (it + 1));
        Lnew(1:2 * p * (it + 1), p + 1:2 * p) = ...
            (I(1:2 * p * (it + 1), 2 * p * it - p + 1:2 * p * it) - ...
             opts.KSM.compute_struct.T(1:2 * p * (it + 1), ...
                                       1:2 * p * it) * ...
             H(1:2 * p * it, p + 1:2 * p)) * ...
             hinv(p + 1:2 * p, p + 1:2 * p);

        opts.KSM.compute_struct.rho = ...
            hinv(1:p, 1:p) \ hinv(1:p, p + 1:2 * p);
        opts.KSM.compute_struct.L = Lnew;

    end

elseif strcmp('RK', opts.KSM.space)

    T = opts.KSM.compute_struct.T;

    % expand T
    if isempty(T)
        AVnew = oper.mul_A(eqn, opts, eqn.type, ...
                           VV(:, p * (it - 1) + 1:p * it), 'N');
        Tnew = VV' * AVnew;
        opts.KSM.compute_struct.T = Tnew;
    else

        AVnew = oper.mul_A(eqn, opts, eqn.type, ...
                           VV(:, p * it + 1:p * (it + 1)), 'N');
        Tnew = VV' * AVnew;

        % compute the new (block) row of T: Vnew^T*A*VV by first computing
        % A^T*Vnew
        % REMARK: in RK, the projection of A is not upper Hessenberg in
        % general
        AVV = oper.mul_A(eqn, opts, eqn.type, VV(:, 1:p * it), 'N');
        Tnew_row = VV(:, p * it + 1:p * (it + 1))' * AVV;
        opts.KSM.compute_struct.T = [[T; Tnew_row], Tnew];
    end
    opts.KSM.compute_struct.AVnew = AVnew;
end

[eqn, opts, oper] = oper.mul_A_post(eqn, opts, oper);
end
