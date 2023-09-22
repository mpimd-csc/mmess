function [out, eqn, opts, oper] = ...
    mess_accumulateK(eqn, opts, oper, out, pc, V1, V2)
% Updates out.Knew and out.DeltaK
%
%  K = E' ZZ' B if eqn.type == 'N'
%  K = E ZZ' C' if eqn.type == 'T'
%
%
% Input:
%  eqn       structure containing equation data
%
%  opts      structure containing parameters for the algorithm
%
%  oper      contains function handles with operations for A and E
%
%  out       contains Knew and DeltaK
%
%  pc        contains shift parameter p
%
%  V1        contains solution of shifted system or Z if pc is empty
%
%  V2        contains solution of shifted system
%
% Output:
%  out       contains Knew and DeltaK
%
%  eqn       structure containing equation data
%
%  opts      structure containing parameters for the algorithm
%
%  oper      contains function handles with operations for A and E

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%% Check input

%% Initialize data
if not(isreal(pc)) && not(nargin == 7)
    mess_err(opts, 'control_data', ...
             'If the shift is complex, V1 and V2 are required');
end

%% preprocess multiplication with E
[eqn, opts, oper] = oper.mul_E_pre(eqn, opts, oper);

%% update K and DeltaK
if isempty(pc)
    if eqn.haveE
        if opts.LDL_T
            % no separate BDF case; the computed K is not the
            % feedback matrix for the ARE resulting in a time
            % step of the BDF method but the feedback matrix of
            % the original DRE; that is why tau and beta do not
            % appear in K, e.g. as factor of B; for residual
            % computations of the ARE tau and beta need to be
            % taken into account.
            if eqn.type == 'T'
                out.Knew = oper.mul_E(eqn, opts, eqn.type, ...
                                      V1 * (out.D * (V1' * eqn.B)), 'N') / ...
                                      eqn.R;
            else
                out.Knew = oper.mul_E(eqn, opts, eqn.type, ...
                                      V1 * (out.D * (eqn.C * V1)'), 'N') / ...
                                      eqn.Q;
            end
        else
            if eqn.type == 'T'
                out.Knew = oper.mul_E(eqn, opts, eqn.type, ...
                                      V1 * (V1' * eqn.B), 'N');
            else
                out.Knew = oper.mul_E(eqn, opts, eqn.type, ...
                                      V1 * (eqn.C * V1)', 'N');
            end
        end
    else
        if opts.LDL_T
            if eqn.type == 'T'
                out.Knew = V1 * (out.D * (V1' * eqn.B)) / eqn.R;
            else
                out.Knew = V1 * (out.D * (eqn.C * V1)') / eqn.Q;
            end
        else
            if eqn.type == 'T'
                out.Knew = V1 * (V1' * eqn.B);
            else
                out.Knew = V1 * (eqn.C * V1)';
            end
        end
    end
else
    if opts.adi.accumulateK || opts.adi.accumulateDeltaK
        if isreal(pc)
            if eqn.haveE
                if eqn.type == 'T'
                    if opts.LDL_T
                        % no separate BDF case; the computed K is not the
                        % feedback matrix for the ARE resulting in a time
                        % step of the BDF method but the feedback matrix of
                        % the original DRE; that is why tau and beta do not
                        % appear in K, e.g. as factor of B; for residual
                        % computations of the ARE tau and beta need to be
                        % taken into account.
                        K_update = ...
                            oper.mul_E(eqn, opts, eqn.type, V1, 'N') * ...
                            ((2 * (-pc) * eqn.T) * (V1' * eqn.B));
                        K_update = K_update / eqn.R;
                    else
                        K_update = ...
                            oper.mul_E(eqn, opts, eqn.type, V1, 'N') * ...
                            ((2 * (-pc)) * (V1' * eqn.B));
                    end
                else
                    if opts.LDL_T
                        K_update = ...
                            oper.mul_E(eqn, opts, eqn.type, V1, 'N') * ...
                            ((2 * (-pc) * eqn.T) * (eqn.C * V1)');
                        K_update = K_update / eqn.Q;
                    else
                        K_update = ...
                            oper.mul_E(eqn, opts, eqn.type, V1, 'N') * ...
                            ((2 * (-pc)) * (eqn.C * V1)');
                    end
                end
            else
                if eqn.type == 'T'
                    if opts.LDL_T
                        K_update = V1 * ((2 * (-pc) * eqn.T) * (V1' * eqn.B));
                        K_update = K_update / eqn.R;
                    else
                        K_update = V1 * ((2 * (-pc)) * (V1' * eqn.B));
                    end
                else
                    if opts.LDL_T
                        K_update = V1 * ((2 * (-pc) * eqn.T) * (eqn.C * V1)');
                        K_update = K_update / eqn.Q;
                    else
                        K_update = V1 * ((2 * (-pc)) * (eqn.C * V1)');
                    end
                end
            end
        else
            if eqn.haveE
                if eqn.type == 'T'
                    if opts.LDL_T
                        K_update = ...
                            oper.mul_E(eqn, opts, eqn.type, V1, 'N') * ...
                            (eqn.T * (V1' * eqn.B)) + ...
                            (oper.mul_E(eqn, opts, eqn.type, V2, 'N')) * ...
                            (eqn.T * (V2' * eqn.B));
                        K_update = K_update / eqn.R;
                    else
                        K_update = ...
                            oper.mul_E(eqn, opts, eqn.type, V1, 'N') * ...
                            (V1' * eqn.B) + ...
                            (oper.mul_E(eqn, opts, eqn.type, V2, 'N')) * ...
                            (V2' * eqn.B);
                    end
                else
                    if opts.LDL_T
                        K_update = ...
                            oper.mul_E(eqn, opts, eqn.type, V1, 'N') * ...
                            (eqn.T * (eqn.C * V1)') + ...
                            (oper.mul_E(eqn, opts, eqn.type, V2, 'N')) * ...
                            (eqn.T * (eqn.C * V2)');
                        K_update = K_update / eqn.Q;
                    else
                        K_update = ...
                            oper.mul_E(eqn, opts, eqn.type, V1, 'N') * ...
                            (eqn.C * V1)' + ...
                            (oper.mul_E(eqn, opts, eqn.type, V2, 'N')) * ...
                            (eqn.C * V2)';
                    end
                end
            else
                if eqn.type == 'T'
                    if opts.LDL_T
                        K_update = V1 * (eqn.T * (V1' * eqn.B)) + ...
                                   V2 * (eqn.T * (V2' * eqn.B));
                        K_update = K_update / eqn.R;
                    else
                        K_update = V1 * (V1' * eqn.B) + V2 * (V2' * eqn.B);
                    end
                else
                    if opts.LDL_T
                        K_update = V1 * (eqn.T * (eqn.C * V1)') + ...
                                   V2 * (eqn.T * (eqn.C * V2)');
                        K_update = K_update / eqn.Q;
                    else
                        K_update = V1 * (eqn.C * V1)' + V2 * (eqn.C * V2)';
                    end
                end
            end
        end
        if opts.adi.accumulateK
            out.Knew = out.Knew + K_update;
        end
        if opts.adi.accumulateDeltaK
            out.DeltaK = out.DeltaK + K_update;
        end
    end
end
%% postprocess multiplication with E
[eqn, opts, oper] = oper.mul_E_post(eqn, opts, oper);
end
