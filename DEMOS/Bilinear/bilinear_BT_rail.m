function bilinear_BT_rail(refinements)
% Demonstrates the basics of truncated Gramian based balanced truncation
% for the bilinear formulation of the rail model.
%
% Call: bilinear_BT_rail(refinements)
%
% Input:
%   refinement   The grid refinement level determining the problem size.
%                See `mess_get_bilinear_rail` for details.  The default
%                uses the size n = 5177 model.
%                (optional, default: 3)

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%% Octave version check
% We need some specific Octave as it seems to ignore mass matrices in old
% versions and many ode functions.
if exist('OCTAVE_VERSION', 'builtin') && ...
        (not(exist('verLessThan', 'file')) || verLessThan('octave', '5.2'))

    disp('This demo needs at least octave 5.2 with ode15s support.');

    return
end

%% Set the input argument, if it is not given.
if nargin < 1
    refinements=3;
end

%% Fetch model data
eqn = mess_get_bilinear_rail(refinements);
%% Set process parameters
% verbosity of the Lyapunov-plus-positive solver
opts.blyap.info = 2;

% Set tolerances and iteration limit for the Lyapunov-plus-positive solver
opts.blyap.res_tol  =  1e-6;
opts.blyap.rel_diff_tol  =  1e-6 ;
opts.blyap.maxiter = 10;

% Set options for the ADI-based inner Lyapunov solver
% and related shift computation
opts.shifts.num_Ritz = 50;
opts.shifts.num_hRitz = 25;
opts.shifts.method = 'projection';
opts.shifts.num_desired = 6;

opts.adi.maxiter = 100;
opts.adi.res_tol = 1e-12;
opts.adi.rel_diff_tol = 0;
opts.adi.info = 2;
opts.adi.norm = 'fro';
opts.norm = 'fro';

% Set mess_res2_norm options
opts.resopts.res.maxiter = 10;
opts.resopts.res.tol = 1e-6;
opts.resopts.res.orth = 0;

opts.fopts.LDL_T = 0;
opts.fopts.norm = 'fro';

opts.srm.tol = 1e-10;
opts.srm.info = 1;
%% choose USFS set
oper = operatormanager('default');

%% Perform system reduction
[ROM, ~, eqn, opts, oper]= mess_balanced_truncation_bilinear(eqn, opts ,oper);
ROM.haveE = 0;
% only needed for function bilinear systems to run for both FOM and ROM:
ROM.A_=ROM.A;
ROM.N_=ROM.N;

%% Evaluation of reduction results
% Set ode45 parameters
t0 = 0;
tf = 4500;
tvals = t0:1:tf;

% Start ode45 for original and reduced order systems
u = ones(size(eqn.B,2),1);
x0 = zeros(size(eqn.A_,1),1);
options=odeset(...
    'RelTol', 1e-6, ...
    'Mass', eqn.E_, ...
    'MStateDependence', 'none',...
    'MassSingular', 'no');

[~, x] = ...
    ode15s(@(t,x)bilinear_system(t,x, u, eqn, opts, oper), tvals, x0, options);

options=odeset('RelTol', 1e-6);
u = ones(size(ROM.B,2),1);
x0 = zeros(1,size(ROM.A,1));
[~,x_r] = ...
    ode15s(@(t,x)bilinear_system(t, x, u, ROM, opts, oper),  tvals, x0, options);

% Calculate and plot system outputs
y   = eqn.C*x';
y_r = ROM.C*x_r';

y = y';
y_r = y_r';

leg_comp = {'FOM out 1', 'FOM out 2', 'FOM out 3', 'FOM out 4', 'FOM out 5', ...
    'FOM out 6', 'ROM out 1', 'ROM out 2', 'ROM out 3', 'ROM out 4',...
    'ROM out 5', 'ROM out 6'}; 
leg_err = {'out 1', 'out 2', 'out 3', 'out 4', 'out 5', 'out 6'};

figure()
plot(tvals, y, '-', 'LineWidth', 3)
hold on
plot(tvals, y_r, '--', 'LineWidth', 3)
xlabel('time [s]')
ylabel('magnitude')
legend(leg_comp, 'Location','northeastoutside')
title('FOM (solid) versus ROM (dashed) outputs')
hold off
% Evaluate absolute error
absErr = abs(y-y_r);
% and corresponding relative error
relErr = abs(absErr ./ y);

% plot relative error
figure()
semilogy(tvals, relErr, 'LineWidth', 3);
xlabel('time [s]')
ylabel('magnitude')
legend(leg_err, 'Location','northeastoutside')
title('pointwise relative output errors')
% Check the relative error
for i = 60:67
    test_condition = relErr(i) < 1e-2;
    if test_condition == 0
        error('limit for relative Error exceeded')
    end
end

fprintf('Everything looks good!\n')

end

%% helper function for the ode integrator
function f = bilinear_system(~, x ,u ,eqn, opts, oper)

    f = eqn.A_*x + eqn.B*u;

    [eqn, opts, oper] = oper.mul_N_pre(eqn, opts, oper);
    numberOf_N_matrices = length(eqn.N_);

    for currentN_k = 1:numberOf_N_matrices
        f = f + oper.mul_N(eqn, opts, 'N', x, 'N', currentN_k)*u(currentN_k);
    end

    [~, ~, ~] = oper.mul_N_post(eqn, opts, oper);

end