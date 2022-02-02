function [rw,  Hp, Hm, Vp, Vm, eqn, opts, oper] = mess_get_ritz_vals(eqn,opts,oper)
% [rw,  Hp, Hm, Vp, Vm, eqn, opts, oper] = mess_get_ritz_vals(eqn,opts,oper)
% 
% Computes a number of Ritz and harmonic Ritz values for the operator
% defined by oper and eqn. The actual numbers are determined by 
%
% opts.shifts.num_desired         total number of values
% opts.shifts.num_Ritz            number of Ritz of values
% opts.shifts.num_hRitz           number of harmonic Ritz values
%
% The initial vector for the Arnoldi processes with "A" or "A^-1" can be
% set in opts.shifts.b0 when unset a vector of all ones will be used.
%

%
% This file is part of the M-M.E.S.S. project 
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%


%% check data
if not(isfield(opts,'shifts')) || not(isstruct(opts.shifts))
    warning('MESS:control_data',['shift parameter control structure missing.', ...
        'Switching to default num_desired = 25, num_Ritz = 50, num_hRitz = 25.']);
    opts.shifts.num_desired = 25;
    opts.shifts.num_Ritz = 50;
    opts.shifts.num_hRitz = 25;
else
    if not(isfield(opts.shifts,'num_desired'))||not(isnumeric(opts.shifts.num_desired))
        warning('MESS:control_data',...
            ['Missing or Corrupted opts.shifts.num_desired field.', ...
            'Switching to default: 25']);
        opts.shifts.num_desired = 25;
    end
    if not(isfield(opts.shifts,'num_Ritz'))||not(isnumeric(opts.shifts.num_Ritz))
        warning('MESS:control_data',...
            ['Missing or Corrupted opts.shifts.num_Ritz field.', ...
            'Switching to default: 50']);
        opts.shifts.num_Ritz = 50;
    end
    if not(isfield(opts.shifts,'num_hRitz'))||not(isnumeric(opts.shifts.num_hRitz))
        warning('MESS:control_data',...
            ['Missing or Corrupted opts.shifts.num_hRitz field.', ...
            'Switching to default: 25']);
        opts.shifts.num_hRitz = 25;
    end
end
if not(isfield(eqn, 'haveE')), eqn.haveE = 0; end
[result, eqn, opts, oper] = oper.init(eqn, opts, oper, 'A','E');
if not(result)
    error('MESS:control_data', 'system data is not completely defined or corrupted');
end

n = oper.size(eqn, opts);

if opts.shifts.num_Ritz >= n, error('num_Ritz must be smaller than n!'); end
if opts.shifts.num_hRitz >= n, error('num_hRitz must be smaller than n!'); end
if (2 * (opts.shifts.num_desired) >= opts.shifts.num_Ritz + opts.shifts.num_hRitz), ...
        error('2*num_desired must be smaller than num_Ritz+num_hRitz!'); end

if (not(isfield(opts.shifts, 'b0')) || isempty(opts.shifts.b0))
    opts.shifts.b0 = ones(n, 1);
end


%% initialize data
opts.shifts.b0 = (1 / norm(opts.shifts.b0)) * opts.shifts.b0;
rw = [];
Hp = [];
Vp = [];
Hm = [];
Vm = [];

%% estimate suboptimal ADI shift parameters
if opts.shifts.num_Ritz > 0
    [Hp, Vp] = mess_arn(eqn, opts, oper, 'N');
    rwp = eig(Hp(1:opts.shifts.num_Ritz, 1:opts.shifts.num_Ritz));      % =: R_+
    rw = [rw; rwp];
end

if opts.shifts.num_hRitz > 0
    [Hm, Vm] = mess_arn(eqn, opts, oper, 'I');
    rwm = ones(opts.shifts.num_hRitz, 1)./eig(Hm(1:opts.shifts.num_hRitz, ...
        1:opts.shifts.num_hRitz));                                  % =: 1 / R_-
    rw = [rw; rwm];                                                       % =: R
end
if any(real(rw) >= zeros(size(rw)))
    warning('MESS:antistable_ritz',...
        ['Non-stable Ritz values were detected.\n',...
        'These will be removed from the set in further computations.']);
    rw  = rw(real(rw)<0);
end
end
