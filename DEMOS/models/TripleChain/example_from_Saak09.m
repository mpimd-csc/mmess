% Simple wrapper script to generate the version of the triple chain
% oscillator model from
% [1] N.Truhar and K.Veselic
%     An efficient method for estimating the optimal dampers' viscosity for
%     linear vibrating systems using Lyapunov equations
%     SIAM J. Matrix Anal. Appl. vol.31 no.1 pp 18-39
%
% with the parametrization used in
%
% [2] J. Saak, Efficient numerical solution of large scale algebraic matrix equations
%     in PDE control and model order reduction, Dissertation,
%     Technische Universität Chemnitz, Chemnitz, Germany (Jul. 2009).
%     URL http://nbn-resolving.de/urn:nbn:de:bsz:ch1-200901642

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%%
n1=150;
alpha=0.01;
beta=alpha;
v=5e0;

%%
[M,E,K]=triplechain_MSD(n1,alpha,beta,v);
B=ones(3*n1+1,1);
Cp=B';
Cv=zeros(size(Cp));


%%
nsample=200;
w=logspace(-4,2,nsample);

tro=zeros(1,nsample);
fprintf(['Computing TFMs of original systems and ' ...
  'MOR errors\n'])

%%
for k=1:nsample
  fprintf('\r Step %3d / %3d',k,nsample)
    Go = (Cp + 1i*w(k)*Cv)/(-w(k)*w(k)*M + 1i*w(k)*D + K) * B;
    tro(k) = max(svds(Go));
end
fprintf('\n\n');
figure(1)
loglog(w, tro, 'LineWidth', 3);
xlabel('\omega');
ylabel('\sigma_{max}(G(j\omega))');
title('Transfer functions of original systems');
