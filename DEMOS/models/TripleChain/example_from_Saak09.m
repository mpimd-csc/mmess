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
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, see <http://www.gnu.org/licenses/>.
%
% Copyright (C) Jens Saak, Martin Koehler, Peter Benner and others 
%               2009-2020
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
%Cv=B';

%%
nsample=200;
%w=logspace(-4,0,nsample);
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
loglog(w, tro)
xlabel('\omega')
ylabel('\sigma_{max}(G(j\omega))')
title('Transfer functions of original systems')
