function [Er, Ar, Br, Cr] = IRKA_rail(k, r, istest)
% Computes a locally H2-optimal reduced order model of order r for
% the selective cooling of Steel profiles application described in
% [1,2,3] via the tangential IRKA method.
%
% Usage: IRKA_rail(k,r,istest)
%
% Inputs:
%
% k           refinement level of the model to use
%             (0 - 5, i.e. 109 - 79841 Dofs)
%             (optional, defaults to 2, i.e. 1357 Dofs)
%
% r           desired dimension of the reduced order model
%             (optional, defaults to 10)
%
% istest      flag to determine whether this demo runs as a CI test or
%             interactive demo
%             (optional, defaults to 0, i.e. interactive demo)
%
% References:
% [1] J. Saak, Effiziente numerische Lösung eines
%     Optimalsteuerungsproblems für die Abkühlung von Stahlprofilen,
%     Diplomarbeit, Fachbereich 3/Mathematik und Informatik, Universität
%     Bremen, D-28334 Bremen (Sep. 2003).
%     https://doi.org/10.5281/zenodo.1187040
%
% [2] P. Benner, J. Saak, A semi-discretized heat transfer model for
%     optimal cooling of steel profiles, in: P. Benner, V. Mehrmann, D.
%     Sorensen (Eds.), Dimension Reduction of Large-Scale Systems, Vol. 45
%     of Lecture Notes in Computational Science and Engineering, Springer-Verlag, Berlin/Heidelberg,
%     Germany, 2005, pp. 353–356. https://doi.org/10.1007/3-540-27909-1_19
%
% [3] J. Saak, Efficient numerical solution of large scale algebraic matrix
%     equations in PDE control and model order reduction, Dissertation,
%     Technische Universität Chemnitz, Chemnitz, Germany (Jul. 2009).
%     URL http://nbn-resolving.de/urn:nbn:de:bsz:ch1-200901642
%
% [4] S. Gugercin, A. C. Antoulas, C. Beattie, H2 model reduction
%     for large-scale linear dynamical systems, SIAM J. Matrix
%     Anal. Appl. 30 (2) (2008) pp. 609–638.
%     https://doi.org/10.1137/060666123

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%% Load matrix data
if nargin < 1
    k = 2;
end
eqn = mess_get_linear_rail(k);

%% register corresponding usfs
opts = struct();
[oper, opts] = operatormanager(opts, 'default');

%% collect IRKA parameters
if nargin < 2

    opts.irka.r = 20;

else

    opts.irka.r = r;

end

if nargin < 3
    istest = false;
end

opts.irka.maxiter   = 100;
opts.irka.shift_tol = 1e-3;
opts.irka.h2_tol    = 1e-6;

if istest

    opts.irka.info = 1;

else

    opts.irka.info = 3;

end

opts.irka.init = 'logspace';

%% Run IRKA

[Er, Ar, Br, Cr, ~, outinfo] = mess_tangential_irka(eqn, opts, oper);

if istest && isequal(outinfo.term_flag, 'maxiter')
    mess_err(opts, 'convergence', 'IRKA converged unexpectedly slow');
end
%% In case this is a CI test and the sparss class is available (recent MATLAB)
%  try the same again with that.
if istest && exist('sparss', 'class')

    sys = sparss(eqn.A_, eqn.B, eqn.C, [], eqn.E_);

    [Er, Ar, Br, Cr, ~, outinfo] = mess_tangential_irka(sys, opts);

    if istest && isequal(outinfo.term_flag, 'maxiter')

        mess_err(opts, 'convergence', 'IRKA converged unexpectedly slow');

    end
end
