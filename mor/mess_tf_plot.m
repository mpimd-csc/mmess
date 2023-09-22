function [out, eqn, opts, oper] = mess_tf_plot(arg_one, opts, varargin)
% Computation of simple  transfer function plots for descriptor systems with
% and comparison to reduced order models.
%
% Calling sequence:
%
%   [out, eqn, opts, oper] = mess_tf_plot(eqn, opts, oper, ROM)
% or
%
%   [out, eqn, opts, oper] = mess_tf_plot(g, opts, ROM)
%
% INPUTS:
%   eqn       struct contains data for equations
%
%   g         presampled transfer function of the original system
%             fitting the parameters in opts.tf_plot
%
%   opts      struct contains parameters for the algorithm
%             (mandatory with substructure opts.tf_plot)
%
%   oper      struct contains function handles for operation with A and E
%
%   ROM       structure containing reduced order model matrices
%             either E, A, B, C, D,
%                 or M, E, K, B, Cv, Cp, D
%             where in the first case E and in both cases D are optional.
%             (optional when eqn is present; mandatory otherwise)
%
% The relevant substructure of opts is opts.tf_plot with members:
%
% type         type of plot (optional, defaults to 'sigma')
%                'sigma'  sigma magnitude plot
%                'Fro'    Frobenius norm plot
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
%  w           vector of frequency samples in rad/s
%              (mandatory for input g, ignored for input eqn)
%
%  Hz          Indicates to use Hertz on the frequency-axis, when info == 2.
%              Only used for plotting, Output frequencies (out.w) will
%              remain given in rad/s.
%              (optional, default: false)
%
%  db          Indicates to use decibels on the magnitude-axis, when
%              info == 2.
%              Only scales presentation in the plot, not the vectors in
%              out below.
%              (optional, default: false)
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
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%% Check and assign inputs
narginchk(2, 4);

if not(isa(opts, 'struct'))
    mess_err(opts, 'illegal_input', ...
             'second input must be am options structure');
end

if nargin == 2
    if not(isa(arg_one, 'numeric'))
        mess_err(opts, 'illegal_input', ...
                 ['First input must be numeric (3d array of doubles)', ...
                  ' in the two input case']);
    else
        g = arg_one;
        ROM = [];
    end
end

if nargin > 2
    if not(isa(arg_one, 'struct')) && not(isa(arg_one, 'numeric'))
        mess_err(opts, 'illegal_input', ...
                 ['input 1 must be numeric (3d array of doubles)', ...
                  ' or equation structure']);
    end

    if not(isa(varargin{1}, 'struct')) || ...
       (nargin == 4 && not(isa(varargin{2}, 'struct')))
        mess_err(opts, 'illegal_input', ...
                 ['Either all inputs are structures, ', ...
                  'or the first input is numeric ', ...
                  'and the rest are structures']);
    end

    if isa(arg_one, 'numeric')
        g = arg_one;
        ROM = varargin{1};
        eqn = [];
        oper = [];
    else

        g = [];
        eqn = arg_one;
        oper = varargin{1};

        if nargin == 4
            ROM = varargin{2};
        else
            ROM = [];
        end
    end
end

%% check field opts.tf_plot
if not(isfield(opts, 'tf_plot')) || not(isstruct(opts.tf_plot))
    mess_err(opts, 'control_data', ...
             'No tf_plot control data found in options structure.');
end

if not(isfield(opts.tf_plot, 'type')) || not(ischar(opts.tf_plot.type))
    mess_warn(opts, 'control_data', ...
              ['Missing or invalid type selector in options ', ...
               'structure. Falling back to sigma magnitude plot.\n']);
    opts.tf_plot.type = 'sigma';
end

%% select plot type and setup axis labels
switch lower(opts.tf_plot.type)
    case 'sigma'
        fun = @(TF) max(svd(TF));
        ystr_err = '\sigma_{max}(G(j\omega) - G_r(j\omega))';
        ystr_fun = '\sigma_{max}(G(j\omega))';

    case {'fro', 'frobenius'}
        fun = @(TF) norm(TF, 'fro');
        ystr_err = '||G(j\omega) - G_r(j\omega)||_{F}';
        ystr_fun = '||G(j\omega)||_{F}';

    otherwise
        mess_err(opts, 'control_data', ...
                 'Unknown plot type requested.');

end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check info parameter for output verbosity and desired axis units
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if not(isfield(opts.tf_plot, 'info'))
    opts.tf_plot.info = 2;
else
    if not(isnumeric(opts.tf_plot.info)) && ...
            not(islogical(opts.tf_plot.info))
        mess_err(opts, 'info', ...
                 'opts.tf_plot.info parameter must be logical or numeric.');
    end
end

if not(isfield(opts.tf_plot, 'Hz'))
    opts.tf_plot.Hz = false;
end
if opts.tf_plot.Hz
    xstr = '\omega [Hz]';
else
    xstr = '\omega [rad/s]';
end

if not(isfield(opts.tf_plot, 'db'))
    opts.tf_plot.db = false;
end
if opts.tf_plot.db
    ystr_fun = [ystr_fun, ' [db]'];
    ystr_err = [ystr_err, ' [db]'];
end
%%
% check if required frequency vector is given together with pre-samples in g
if isempty(g)

    if not(isfield(opts.tf_plot, 'w')) && ...
       not(all(isfield(opts.tf_plot, {'fmin', 'fmax'})))

        mess_err(opts, 'control_data', ...
                 ['tf_plot control data does not contain', ...
                  ' frequency range bounds.']);
    end

    if not(isfield(opts.tf_plot, 'nsample'))

        opts.tf_plot.nsample = 100;
    end

elseif not(isfield(opts.tf_plot, 'w'))

    mess_err(opts, 'control_data', ...
             ['tf_plot control data must contain frequency vector w', ...
              ' when presampled original is given']);
end

% if no frequency vector was passed in let's generate one
if isfield(opts.tf_plot, 'w')
    w = opts.tf_plot.w;
    opts.tf_plot.nsample = length(w);
else
    if (floor(opts.tf_plot.fmin) == opts.tf_plot.fmin) && ...
       (floor(opts.tf_plot.fmax) == opts.tf_plot.fmax)
        w = logspace(opts.tf_plot.fmin, ...
                     opts.tf_plot.fmax, ...
                     opts.tf_plot.nsample);
    else
        w = logspace(log10(opts.tf_plot.fmin), ...
                     log10(opts.tf_plot.fmax), ...
                     opts.tf_plot.nsample);
    end
end

%% preallocation
out.tr1 = zeros(1, opts.tf_plot.nsample);

if isempty(g)

    m = size(eqn.C, 1);
    p = size(eqn.B, 2);

    out.g1 = zeros(m, p, opts.tf_plot.nsample);
else

    [m, p] = size(g(:, :, 1));

    out.g1 = g;
end

if not(isempty(ROM))
    out.tr2 = out.tr1;
    out.err = out.tr1;
    out.relerr = out.tr1;

    out.g2 = zeros(m, p, opts.tf_plot.nsample);
end

if opts.tf_plot.info
    mess_fprintf(opts, ...
                 ['Computing TFMs of original and reduced order systems ' ...
                  'and MOR errors\n']);
end

%% preprocess shifted solver if eqn is given
if isempty(g)
    [result, eqn, opts, oper] = oper.init(eqn, opts, oper, 'A', 'E');

    if not(result)
        mess_err(opts, 'control_data', ...
                 'system data is not completely defined or corrupted');
    end

    [eqn, opts, oper] = oper.sol_ApE_pre(eqn, opts, oper);
end

%% make sure we have an E in the first-order ROM
if not(isempty(ROM))

    if (not(isfield(ROM, 'A')) && any(not(isfield(ROM, {'K', 'M'})))) || ...
            not(isfield(ROM, 'B')) || ...
            (not(isfield(ROM, 'C')) && ...
             all(not(isfield(ROM, {'Cv', 'Cp'}))))

        mess_err(opts, 'illegal_input', ...
                 'Found incomplete ROM structure!');
    end

    if not(isfield(ROM, 'E'))
        if not(isfield(ROM, 'A'))
            ROM.E = [];
        else
            ROM.E = eye(size(ROM.A));
        end
    end

    if isfield(ROM, 'Cp') && not(isfield(ROM, 'Cv'))
        ROM.Cv = [];
    end
    if isfield(ROM, 'Cv') && not(isfield(ROM, 'Cp'))
        ROM.Cp = [];
    end
end

%% perform the actual sampling
for k = 1:opts.tf_plot.nsample

    if opts.tf_plot.info && not(mod(k, opts.tf_plot.nsample / 10))
        mess_fprintf(opts, '\r Step %3d / %3d', k, opts.tf_plot.nsample);
    end

    if isempty(g) % sample original model only if it was not given

        out.g1(:, :, k) = full(eqn.C * ...
                               oper.sol_ApE(eqn, opts, ...
                                            'N', -1i * w(k), ...
                                            'N', -eqn.B, 'N'));

        if isfield(eqn, 'D') && not(isempty(eqn.D))

            out.g1(:, :, k) = out.g1(:, :, k) + eqn.D;
        end
    end

    if not(isempty(ROM)) % sample reduced model only if it was given

        if isfield(ROM, 'A')
            out.g2(:, :, k) = ROM.C * ((1i * w(k) * ROM.E - ROM.A) \ ROM.B);
        else
            out.g2(:, :, k) = (ROM.Cp + 1i * w(k) * ROM.Cv) * ...
                ((-w(k) * (w(k) * ROM.M - 1i * ROM.E) + ROM.K) \ ROM.B);
        end

        if isfield(ROM, 'D') && not(isempty(ROM.D))

            out.g2(:, :, k) = out.g2(:, :, k) + ROM.D;
        end

    end

    out.tr1(k) = max(svd(out.g1(:, :, k)));

    if not(isempty(ROM))
        out.err(k) = fun(out.g1(:, :, k) - out.g2(:, :, k));
        out.tr2(k) = fun(out.g2(:, :, k));
        out.relerr(k) = out.err(k) / out.tr1(k);
    end
end

out.w = w;

if opts.tf_plot.info
    mess_fprintf(opts, '\n\n');
end

%% postprocess shifted solver if eqn is given
if isempty(g)
    [eqn, opts, oper] = oper.sol_ApE_post(eqn, opts, oper);
end

%% Finally,  the plots (if desired)
if opts.tf_plot.Hz
    % we want to rescale to Hz
    rescale_to_unit = @(w) w / (2 * pi);
else
    % do nothing
    rescale_to_unit = @(w) w;
end

if isnumeric(opts.tf_plot.info) && opts.tf_plot.info > 1
    if not(isempty(ROM))

        figure();

        subplot(2, 1, 1);
        if opts.tf_plot.db
            semilogx(rescale_to_unit(out.w), ...
                     20 * log10(squeeze(out.err)), ...
                     'LineWidth', 3);
        else
            loglog(rescale_to_unit(out.w), out.err, 'LineWidth', 3);
        end
        title('absolute model reduction error');
        xlabel(xstr);
        ylabel(ystr_err);
        axis tight;

        subplot(2, 1, 2);
        if opts.tf_plot.db
            semilogx(rescale_to_unit(out.w), ...
                     20 * log10(squeeze(out.relerr)), ...
                     'LineWidth', 3);
        else
            loglog(rescale_to_unit(out.w), out.relerr, 'LineWidth', 3);
        end
        title('relative model reduction error');
        xlabel(xstr);
        ylabel([ystr_err, '/ ', ystr_fun]);
        axis tight;
    end

    figure();
    if opts.tf_plot.db
        semilogx(rescale_to_unit(out.w), ...
                 20 * log10(squeeze(out.tr1)), ...
                 'LineWidth', 3);
    else
        loglog(rescale_to_unit(out.w), out.tr1, 'LineWidth', 3);
    end

    if not(isempty(ROM))
        hold on;
        if opts.tf_plot.db
            semilogx(rescale_to_unit(out.w), ...
                     20 * log10(squeeze(out.tr2)), ...
                     'LineWidth', 3);
        else
            loglog(rescale_to_unit(out.w), out.tr2, 'r--', 'LineWidth', 3);
        end

        hold off;

        legend({'original system', 'reduced system'});
        title('Transfer functions of original and reduced systems');
    else
        legend({'original system'});
        title('Transfer functions of original system');
    end

    xlabel(xstr);
    ylabel(ystr_fun);
    axis tight;
end
