function [eqn, opts, oper] = KSM_space_expansion(eqn, opts, oper)
% function compute_struct = KSM_space_expansion(eqn, opts, oper, ...
%                                                     compute_struct)
% Function that computes the next basis block for the space
% enlarging the basis V and update the projection of A onto such a
% space
%
% Input and output:
%
%   eqn       structure containing equation data
%
%   opts      structure containing parameters for the algorithm
%
%   oper      contains function handles with operations for A and E
%
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
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

p = size(eqn.W, 2);

VV = opts.KSM.compute_struct.V;

% current number of iteration
if strcmp('EK', opts.KSM.space)
    [eqn, opts, oper] = oper.mul_A_pre(eqn, opts, oper);
    [eqn, opts, oper] = oper.sol_A_pre(eqn, opts, oper);
    opts.KSM.compute_struct.it = size(VV, 2) / (2 * p);
elseif strcmp('RK', opts.KSM.space)
    opts.KSM.compute_struct.it = size(VV, 2) / p;
end

if isempty(VV)
    if strcmp('EK', opts.KSM.space)
        rhs = eqn.ssW;
        rhs_inv = oper.sol_A(eqn, opts, eqn.type, rhs, 'N');
        [Vnew, beta] = qr([rhs, rhs_inv], 0);
        opts.KSM.compute_struct.H = [];
        opts.KSM.compute_struct.T = [];
        opts.KSM.compute_struct.L = [];
        opts.KSM.compute_struct.beta = beta;
        opts.KSM.compute_struct.V = Vnew;

        % Compute also the projection of the quadratic term core factor if
        % we are solving a care
        if strcmp(opts.KSM.type_eqn, 'CARE')
            switch eqn.type
                case 'N'
                    opts.KSM.compute_struct.Bm = Vnew' * eqn.ssC';

                case 'T'
                    opts.KSM.compute_struct.Bm = Vnew' * eqn.ssB;
            end
        end
    elseif strcmp('RK', opts.KSM.space)
        rhs = eqn.ssW;

        % orthogonalize the low-rank factor of the rhs to get the initial
        % block
        [Vnew, beta] = qr(rhs, 0);
        opts.KSM.compute_struct.H = [];
        opts.KSM.compute_struct.T = [];
        opts.KSM.compute_struct.L = [];
        opts.KSM.compute_struct.beta = beta;
        opts.KSM.compute_struct.V = Vnew;
        opts.KSM.compute_struct.it = 1;

        % update the projection of A onto the current subspace
        [eqn, opts, oper] = KSM_getT(eqn, opts, oper);

        % Compute the first shift
        opts.KSM.compute_struct.shifts = opts.KSM.init_shifts(2);
        [eqn, opts, oper] = KSM_compute_shifts(eqn, opts, oper);

        % Compute also the projection of the quadratic terms core factor if
        % we are solving a care
        if strcmp(opts.KSM.type_eqn, 'CARE')
            switch eqn.type
                case 'N'
                    opts.KSM.compute_struct.Bm = Vnew' * eqn.ssC';

                case 'T'
                    opts.KSM.compute_struct.Bm = Vnew' * eqn.ssB;
            end
        end

    end

else

    if strcmp('EK', opts.KSM.space)

        % take the last basis block
        last_basis_block = VV(:, end - 2 * p + 1:end);
        % multiply the first p columns by A
        V1 = oper.mul_A(eqn, opts, eqn.type, ...
                        last_basis_block(:, 1:p), 'N');
        % multiply the last p columns by A^{-1}
        V2 = oper.sol_A(eqn, opts, eqn.type, ...
                        last_basis_block(:, p + 1:2 * p), 'N');

        opts.KSM.compute_struct.Vnew = [V1, V2];

        % Block modified Gram-Schmidt
        [eqn, opts, oper] = KSM_mgs(eqn, opts, oper);

        % update the projection of A onto the current subspace
        [eqn, opts, oper] = KSM_getT(eqn, opts, oper);

    elseif strcmp('RK', opts.KSM.space)

        % we want a real basis!
        opts.KSM.compute_struct.cmplxconjugate_flag = false;

        while not(opts.KSM.compute_struct.cmplxconjugate_flag)
            % take the last basis block
            last_basis_block = VV(:, end - p + 1:end);

            % solve the shifted linear system
            % REMARK: mess_solve_shifted_system.m solves (A+p*E)x=v while
            % in the RKSM framework we usually work with (A-p*E)x=v.
            % Therefore, we pass -p as input
            [Vnew, eqn, opts, oper] = ...
                mess_solve_shifted_system(eqn, opts, oper, ...
                                          -opts.KSM.compute_struct.new_shift, ...
                                          last_basis_block);

            opts.KSM.compute_struct.basis_cmplxflag = false;
            if any(imag(Vnew(:))) && ...
                    not(opts.KSM.compute_struct.basis_cmplxflag)
                Vnew = real(Vnew);
                opts.KSM.compute_struct.basis_cmplxflag = true;
            elseif any(imag(Vnew(:))) && ...
                    opts.KSM.compute_struct.basis_cmplxflag
                Vnew = imag(Vnew);
            end

            opts.KSM.compute_struct.Vnew = Vnew;

            % Block modified Gram-Schmidt
            [eqn, opts, oper] = KSM_mgs(eqn, opts, oper);

            if opts.KSM.compute_struct.shift_cmplxflag
                % use the complex conjugate of the current shift as new
                % shift
                opts.KSM.compute_struct.new_shift = ...
                    conj(opts.KSM.compute_struct.new_shift);
                opts.KSM.compute_struct.shifts = ...
                    [opts.KSM.compute_struct.shifts, ...
                     opts.KSM.compute_struct.new_shift];
                opts.KSM.compute_struct.shift_cmplxflag = false;

                % update the projection of A onto the current subspace
                [eqn, opts, oper] = KSM_getT(eqn, opts, oper);

                % update the number of iterations
                VV = opts.KSM.compute_struct.V;
                opts.KSM.compute_struct.it = size(VV, 2) / p;

            else
                opts.KSM.compute_struct.cmplxconjugate_flag = true;
            end
        end

        % Compute next shift
        [eqn, opts, oper] = KSM_compute_shifts(eqn, opts, oper);

        % update the projection of A onto the current subspace
        [eqn, opts, oper] = KSM_getT(eqn, opts, oper);

    end
end

if strcmp('EK', opts.KSM.space)
    [eqn, opts, oper] = oper.mul_A_post(eqn, opts, oper);
    [eqn, opts, oper] = oper.sol_A_post(eqn, opts, oper);
end

end
