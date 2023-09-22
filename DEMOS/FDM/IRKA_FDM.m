function [Er, Ar, Br, Cr] = IRKA_FDM(n0, r, istest)
% IRKA_FDM computes a reduced order model via the IRKA method (see
% e.g. [1]) for a finite difference discretized convection
% diffusion model on the unit square described in [2].
%
% Usage:
%    [Ar, Br, Cr] = IRKA_FDM(n0, r, test)
%
% Inputs
%
% n0          n0^2 gives the dimension of the original model, i.e. n0 is
%             the number of degrees of freedom, i.e. grid points, per
%             spatial direction
%             (optional; defaults to 50)
%
% r           desired dimension of the reduced order model
%             (optional, defaults to 10)
%
% istest      flag to determine whether this demo runs as a CI test or
%             interactive demo
%             (optional, defaults to 0, i.e. interactive demo)
%
% Outputs
%
% Ar, Br, Cr  the reduced order system matrices.
%
% References
% [1] A. C. Antoulas, Approximation of Large-Scale Dynamical Systems, Vol.
%     6 of Adv. Des. Control, SIAM Publications, Philadelphia, PA, 2005.
%     https://doi.org/10.1137/1.9780898718713
%
% [2] T. Penzl, LyaPack Users Guide, Tech. Rep. SFB393/00-33,
%     Sonderforschungsbereich 393 Numerische Simulation auf massiv
%     parallelen Rechnern, TU Chemnitz, 09107 Chemnitz, Germany,
%     available from http://www.tu-chemnitz.de/sfb393/sfb00pr.html (2000).

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%%
if nargin < 1
    n0 = 100;
end

E = speye(n0^2);
A = fdm_2d_matrix(n0, '10*x', '100*y', '0');
B = fdm_2d_vector(n0, '.1<x<=.3');
C = fdm_2d_vector(n0, '.7<x<=.9')';

%%
if nargin < 2
    opts.irka.r = 10;
else
    opts.irka.r = r;
end
if nargin < 3
    istest = false;
end
opts.irka.maxiter = 20;
opts.irka.shift_tol = 1e-2;
opts.irka.h2_tol = 1e-6;
if istest
    opts.irka.info = 1;
else
    opts.irka.info = 2;
end
opts.irka.init = 'logspace';

[Er, Ar, Br, Cr, ~, outinfo] = mess_tangential_irka(E, A, B, C, [], opts);

if istest && isequal(outinfo.term_flag, 'maxiter')
    mess_err(opts, 'convergence', 'IRKA converged unexpectedly slow');
end
