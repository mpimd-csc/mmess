function [err, w, eqn, opts, oper] = mess_sigma_plot(eqn, opts, oper, ROM)
% Computation of simple sigma-plots for descriptor systems with invertible
% E.
%
%   mess_sigma_plot(eqn, opts, oper, ROM)
%
% INPUTS:
%   eqn                 struct contains data for equations
%
%   opts                struct contains parameters for the algorithm
%
%   oper                struct contains function handles for operation 
%                       with A and E
% 
%   ROM                 structure containing reduced order model matrices
%                       either E, A, B, C, D, 
%                       or M, E, K, B, Cv, Cp, D
%  
% The relevant substructure of opts is opts.sigma with members:
%
% fmin, fmax   left and right bounds of the frequency range. They will be
%              interpreted as exponents in the logarithmic range if
%              integers are passed.
%
% nsample      number of transfer function samples to take in the plot
%              (optional, defaults to 100)
%
% info         verbosity control. (optional, defaults to 2)
%                1   only progress is reported
%                2   the actual sigma plots.
%
% Outputs:
%
% w            vector of sampling frequencies 
%
% err          vector sampled error values
%

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
%               2009-2019
%

%% check field opts.sigma
if not(isfield(opts,'sigma')) || not(isstruct(opts.sigma))
    error('MESS:control_data',...
        'No sigma plot control data found in options structure.');
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check info parameter for output verbosity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if not(isfield(opts.sigma,'info'))
    opts.sigma.info=2;
else
    if not(isnumeric(opts.sigma.info))&&not(islogical(opts.sigma.info))
        error('MESS:info',...
            'opts.sigma.info parameter must be logical or numeric.');
    end
end

if not(isfield(opts.sigma,'fmin')) || not(isfield(opts.sigma,'fmax'))
    error('MESS:control_data',...
        'sigma plot control data does not contain frequency rang bounds.');
end

if not(isfield(opts.sigma,'nsample'))
    opts.sigma.nsample = 100;
end

if (floor(opts.sigma.fmin)==opts.sigma.fmin) && (floor(opts.sigma.fmax)==opts.sigma.fmax)
    w=logspace(opts.sigma.fmin,opts.sigma.fmax,opts.sigma.nsample);
else
    w=logspace(log10(opts.sigma.fmin),log10(opts.sigma.fmax),opts.sigma.nsample);     
end

tr1=zeros(1,opts.sigma.nsample); tr2=tr1; err=tr1; relerr=tr1;

fprintf(['Computing TFMs of original and reduced order systems and ' ...
         'MOR errors\n']);
%% preprocess shifted solver

[result, eqn, opts, oper] = oper.init(eqn, opts, oper, 'A','E');
if not(result)
    error('MESS:control_data', 'system data is not completely defined or corrupted');
end

[eqn, opts, oper] = oper.sol_ApE_pre(eqn, opts, oper);

%% make sure we have an E in the first-order ROM
if isfield(ROM,'A')  && not(isfield(ROM,'E'))
    ROM.E = eye(size(ROM.A));
end

%% perform the actual sampling
for k=1:opts.sigma.nsample
  if (opts.sigma.info && not(mod(k,opts.sigma.nsample/10)))
      fprintf('\r Step %3d / %3d',k,opts.sigma.nsample); 
  end
  if isfield(eqn,'D')&& not(isempty(eqn.D))
      g1 = full(eqn.C *oper.sol_ApE(eqn, opts,'N',-1i*w(k),'N',-eqn.B,'N')+eqn.D);
  else
      g1 = full(eqn.C *oper.sol_ApE(eqn, opts,'N',-1i*w(k),'N',-eqn.B,'N'));
  end
  if isfield(ROM,'D') && not(isempty(ROM.D))
      if isfield(ROM,'A')
         g2 = ROM.C * ( (1i*w(k)*ROM.E - ROM.A) \ ROM.B ) + ROM.D;
      else
         g2 = (ROM.Cp + 1i*w(k)*ROM.Cv) * (( -w(k) * ( w(k)*ROM.M - 1i*ROM.E ) + ROM.K) \ ROM.B ) + ROM.D;
      end
  else
      if isfield(ROM,'A')
        g2 = ROM.C * ( (1i*w(k)*ROM.E - ROM.A) \ ROM.B );
      else
        g2 = (ROM.Cp + 1i*w(k)*ROM.Cv) * (( -w(k) * ( w(k)*ROM.M - 1i*ROM.E ) + ROM.K) \ ROM.B );    
      end
  end
  err(k) = max(svd(g1-g2));
  tr1(k) = max(svd(g1));
  tr2(k) = max(svd(g2));
  relerr(k)=err(k)/tr1(k);
end
fprintf('\n\n');

%% postprocess shifted solver
[eqn, opts, oper] = oper.sol_ApE_post(eqn, opts, oper); 

if isnumeric(opts.sigma.info) && opts.sigma.info > 1
    figure;
    subplot(2,1,1);
    loglog(w, err);
    title('absolute model reduction error');
    xlabel('\omega');
    ylabel('\sigma_{max}(G(j\omega) - G_r(j\omega))');
    axis tight;
    subplot(2,1,2);
    loglog(w, relerr);
    title('relative model reduction error');
    xlabel('\omega');
    ylabel(['\sigma_{max}(G(j\omega) - G_r(j\omega)) / \' ...
        'sigma_{max}(G(j\omega))']);
    axis tight;
    
    figure;
    loglog(w, tr1);
    hold on;
    loglog(w, tr2, 'r--');
    legend({'original system','reduced system'});
    xlabel('\omega');
    ylabel('\sigma_{max}(G(j\omega))');
    title('Transfer functions of original and reduced systems');
    axis tight;
    hold off;
end
