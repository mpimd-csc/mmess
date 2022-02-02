function [out, eqn, opts, oper] = mess_sigma_plot(arg_one, opts, varargin)
% Computation of simple sigma-magnitude-plots for descriptor systems with
% invertible E and comparison to reduced order models.
%
% Calling sequence:
%
%   [out, eqn, opts, oper] = mess_sigma_plot(eqn, opts, oper, ROM)
% or
%   [out, eqn, opts, oper] = mess_sigma_plot(g, opts, ROM)
%
% INPUTS:
%   eqn                 struct contains data for equations
%
%   g                   presampled transfer function of the original system
%                       fitting the parameters in opts.sigma
%
%   opts                struct contains parameters for the algorithm
%                       (mandatory with substructure opts.sigma)
%
%   oper                struct contains function handles for operation
%                       with A and E
%
%   ROM                 structure containing reduced order model matrices
%                       either E, A, B, C, D,
%                       or M, E, K, B, Cv, Cp, D
%                       where in the first case E and in both cases D are
%                       optional.
%                       (optional when eqn is present; mandatory otherwise)
%
% The relevant substructure of opts is opts.sigma with members:
%
% fmin, fmax   left and right bounds of the frequency range. They will be
%              interpreted as exponents in the logarithmic range if
%              integers are passed.
%              (mandatory for input eqn, ignored for input tr (check w
%              below))
%
% nsample      number of transfer function samples to take in the plot
%              (optional, defaults to 100)
%
% info         verbosity control. (optional, defaults to 2)
%                1   only progress is reported
%                2   also generate the actual sigma plots.
%
%  w           vector of frequency samples fitting
%              (mandatory for input tr, ignored for input eqn)
%
% Outputs:
%
% out.w         vector of sampling frequencies
%
% out.err       vector of sampled maximal singular values of the transfer
%               function of the error system (only when ROM was given)
% out.relerr    vector of sampled maximal singular values of the transfer
%               function of the error system (only when ROM was given)
%               relative to the FOM.
% out.tr1       vector of sampled maximal singular values of the transfer
%               function of the FOM
% out.tr2       vector of sampled maximal singular values of the transfer
%               function of the ROM (only when ROM is given)
% out.g1        vector of sampled values of the transfer function of the
%               FOM
% out.g2        vector of sampled values of the transfer function of the
%               ROM (only when ROM is given)


%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%


%% Check and assign inputs
narginchk(2,4);

if not(isa(opts,'struct'))
    error('second input must be am options structure');
end

if nargin == 2
    if not(isa(arg_one,'numeric'))
        error(['input 1 must be numeric (3d array of doubles)',...
               ' in the 2 input case']);
    else
        g = arg_one;
        ROM = [];
    end
end

if nargin > 2
    if (not(isa(arg_one, 'struct')) && not(isa(arg_one, 'numeric')))
        error(['input 1 must be numeric (3d array of doubles)',...
               ' or equation structure']);
    end

    if not(isa(varargin{1}, 'struct')) || ...
       (nargin == 4 && not(isa(varargin{2}, 'struct')))
        error(['Either all inputs are structures, or the first input',...
               'is numeric and the rest are structures'])
    end

    if isa(arg_one, 'numeric')
        g = arg_one;
        ROM = varargin{1};
        eqn = [];
        oper = [];
    else

        g=[];
        eqn = arg_one;
        oper = varargin{1};

        if nargin == 4
            ROM = varargin{2};
        else
            ROM = [];
        end
    end
end

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
    opts.sigma.info = 2;
else
    if not(isnumeric(opts.sigma.info)) && not(islogical(opts.sigma.info))
        error('MESS:info', ...
              'opts.sigma.info parameter must be logical or numeric.');
    end
end

if isempty(g)

    if not(isfield(opts.sigma,'w')) && ...
       not(all(isfield(opts.sigma,{'fmin','fmax'})))

        error('MESS:control_data',...
              ['sigma plot control data does not contain', ...
               ' frequency range bounds.']);
    end

    if not(isfield(opts.sigma,'nsample'))

        opts.sigma.nsample = 100;
    end

elseif not(isfield(opts.sigma,'w'))

    error('MESS:control_data', ...
          ['sigma plot control data must contain frequency vector w', ...
           ' when presampled original is given']);
end

if isfield(opts.sigma,'w')
    w = opts.sigma.w;
    opts.sigma.nsample = length(w);
else
    if (floor(opts.sigma.fmin) == opts.sigma.fmin) && ...
       (floor(opts.sigma.fmax) == opts.sigma.fmax)
        w = logspace(opts.sigma.fmin,opts.sigma.fmax,opts.sigma.nsample);
    else
        w = logspace(log10(opts.sigma.fmin), ...
                     log10(opts.sigma.fmax),opts.sigma.nsample);
    end
end

%% preallocation
out.tr1 = zeros(1,opts.sigma.nsample);

if isempty(g)

    m = size(eqn.C,1);
    p = size(eqn.B,2);

    out.g1 = zeros(m,p,opts.sigma.nsample);
else

    [m,p] = size(g(:,:,1));

    out.g1 = g;
end

if not(isempty(ROM))
    out.tr2 = out.tr1;
    out.err = out.tr1;
    out.relerr = out.tr1;

    out.g2 = zeros(m,p,opts.sigma.nsample);
end

if opts.sigma.info
    fprintf(['Computing TFMs of original and reduced order systems and ' ...
             'MOR errors\n']);
end

%% preprocess shifted solver if eqn is given
if isempty(g)
    [result, eqn, opts, oper] = oper.init(eqn, opts, oper, 'A','E');

    if not(result)
        error('MESS:control_data', ['system data is not completely',...
              ' defined or corrupted']);
    end

    [eqn, opts, oper] = oper.sol_ApE_pre(eqn, opts, oper);
end

%% make sure we have an E in the first-order ROM
if not(isempty(ROM))

    if (not(isfield(ROM,'A')) && any(not(isfield(ROM,{'K','M'})))) || ...
       not(isfield(ROM,'B')) || ...
       (not(isfield(ROM,'C')) && all(not(isfield(ROM,{'Cv','Cp'}))))
         error('Found incomplete ROM structure!')
    end

    if not(isfield(ROM,'E'))
        if not(isfield(ROM,'A'))
            ROM.E = [];
        else
            ROM.E = eye(size(ROM.A));
        end
    end

    if isfield(ROM, 'Cp') && not(isfield(ROM, 'Cv')), ROM.Cv = []; end
    if isfield(ROM, 'Cv') && not(isfield(ROM, 'Cp')), ROM.Cp = []; end
end

%% perform the actual sampling
for k=1:opts.sigma.nsample

    if (opts.sigma.info && not(mod(k,opts.sigma.nsample/10)))
        fprintf('\r Step %3d / %3d',k,opts.sigma.nsample);
    end

    if isempty(g) % sample original model only if it was not given

        out.g1(:,:,k) = full(eqn.C * oper.sol_ApE(eqn, opts, ...
                             'N',-1i*w(k),'N',-eqn.B,'N'));

        if isfield(eqn,'D') && not(isempty(eqn.D))

            out.g1(:,:,k) = out.g1(:,:,k) + eqn.D;
        end
    end

    if not(isempty(ROM)) % sample reduced model only if it was given

        if isfield(ROM,'A')
            out.g2(:,:,k) = ROM.C * ( (1i*w(k)*ROM.E - ROM.A) \ ROM.B );
        else
            out.g2(:,:,k) = (ROM.Cp + 1i*w(k)*ROM.Cv) * ...
                (( -w(k) * ( w(k)*ROM.M - 1i*ROM.E ) + ROM.K) \ ROM.B );
        end

        if isfield(ROM,'D') && not(isempty(ROM.D))

            out.g2(:,:,k) = out.g2(:,:,k) + ROM.D;
        end

    end

    out.tr1(k) = max(svd(out.g1(:,:,k)));

    if not(isempty(ROM))
        out.err(k) = max(svd(out.g1(:,:,k) - out.g2(:,:,k)));
        out.tr2(k) = max(svd(out.g2(:,:,k)));
        out.relerr(k) = out.err(k)/out.tr1(k);
    end
end

out.w = w;

if opts.sigma.info
    fprintf('\n\n');
end

%% postprocess shifted solver if eqn is given
if isempty(g)
    [eqn, opts, oper] = oper.sol_ApE_post(eqn, opts, oper);
end

%% Finally,  the plots (if desired)
if isnumeric(opts.sigma.info) && opts.sigma.info > 1
    if not(isempty(ROM))

        figure();

        subplot(2,1,1);
        loglog(out.w, out.err, 'LineWidth', 3);
        title('absolute model reduction error');
        xlabel('\omega');
        ylabel('\sigma_{max}(G(j\omega) - G_r(j\omega))');
        axis tight;

        subplot(2,1,2);
        loglog(out.w, out.relerr, 'LineWidth', 3);
        title('relative model reduction error');
        xlabel('\omega');
        ylabel(['\sigma_{max}(G(j\omega) - G_r(j\omega)) / \' ...
                'sigma_{max}(G(j\omega))']);
        axis tight;
    end

    figure();

    loglog(out.w, out.tr1, 'LineWidth', 3);

    if not(isempty(ROM))
        hold on;
        loglog(out.w, out.tr2, 'r--', 'LineWidth', 3);
        hold off;

        legend({'original system','reduced system'});
        title('Transfer functions of original and reduced systems');
    else
        legend({'original system'});
        title('Transfer functions of original system');
    end

    xlabel('\omega');
    ylabel('\sigma_{max}(G(j\omega))');
    axis tight;
end

