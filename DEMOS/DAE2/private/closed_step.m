function closed_step(eqn,Ar,Br,Cr,problem,K,Kr,istest)
% Simple validation of the DAE2 MESS closed loop example via a basic
% step response computation
%
%  ope_step(eqn,Ar,Br,Cr,problem, istest)
%
% INPUTS:
% eqn          The original system equation structure
% Ar,Br,Cr     The reduced order system matrices
% K,Kr         full and reduced order feedback matrices
% problem      'NSE' or 'Stokes' switching between the Stokes demo or the
%               linearized Navier-Stokes-Equation
% istest        flag to determine whether this demo runs as a CI test or
%               interactive demo

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

x0=zeros(size(eqn.A_,1),1);
xr0=zeros(size(Ar,1),1);

alpha=1;
tau=1e-2;
tmin=0;
tmax=50;
T=tmin:tau:tmax;
ntau=floor((tmax-tmin)/tau);
taurange=1:floor(ntau/500):ntau;
%%
t_impeuler_closed = tic;
[y,yr] = impeuler_closed(eqn.E_,eqn.A_,eqn.B,eqn.C,eye(size(Ar)),Ar,Br,...
    Cr,K,Kr,tau,tmin,tmax,x0,xr0,alpha);
t_elapsed = toc(t_impeuler_closed);
fprintf(1,'implicit euler took %6.2f seconds \n', t_elapsed);
%%
abserr=abs(y-yr);
relerr=abs(abserr./y);

%%
if istest
    maerr = max(abserr);
    if maerr>=1e-6
        error('MESS:TEST:accuracy',['unexpectedly inaccurate result ' ...
                            'in closed loop simulation. Maximum ' ...
                            'absolute error %e > 1e-6'], maerr);

    end
else
    colors=['y','m','c','r','g','b','k'];

    figure();
    hold on;
    for j=1:size(eqn.C,1)
        plot(T(taurange),y(j,taurange),colors(j),'LineWidth',3);
        plot(T(taurange),yr(j,taurange),strcat(colors(j),'--'),'LineWidth',3);
    end
    xlabel('time');
    ylabel('magnitude of outputs');
    title('step response');
    if strcmp(problem,'NSE')
        legend('out1','out1 red','out2','out2 red','out3','out3 red',...
            'out4','out4 red','out5','out5 red','out6','out6 red',...
            'out7','out7 red','Location','EastOutside');
    else
        legend('out1','out1 red','out2','out2 red','out3','out3 red',...
            'out4','out4 red','out5','out5 red','Location','EastOutside');
    end
    hold off;
    %%
    figure();
    for j=1:size(eqn.C,1)
        semilogy(T(taurange),abserr(j,taurange),colors(j),'LineWidth',3);
        if 1==j, hold on; end
    end
    xlabel('time');
    ylabel('magnitude');
    title('absolute error');
    if strcmp(problem,'NSE')
        legend('out1','out2','out3','out4','out5','out6','out7',...
        'Location','EastOutside');
    else
        legend('out1','out2','out3','out4','out5','Location','EastOutside');
    end
    hold off;

    figure();

    for j=1:size(eqn.C,1)
        semilogy(T(taurange),relerr(j,taurange),colors(j),'LineWidth',3);
        if 1==j, hold on; end
    end
    xlabel('time');
    ylabel('magnitude');
    title('relative error');
    if strcmp(problem,'NSE')
        legend('out1','out2','out3','out4','out5','out6','out7',...
            'Location','EastOutside');
    else
        legend('out1','out2','out3','out4','out5','Location','EastOutside');
    end
    hold off;
end
