function [out, eqn, opts, oper] = mess_sigma_plot(arg_one, opts, varargin)
% Computation of simple sigma-magnitude-plots for descriptor systems and
% comparison to reduced order models.
%
% Backward compatibility wrapper for mess_tf_plot with
%    opts.tf_plot.type == 'sigma'
%
% arguments are as in mess_tf_plot. All settings put into opts.sigma are
% forwarded to opts.tf_plot.
%

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%% For backward compatibility copy all information from
%  opts.sigma to opts.tf_plot
if isfield(opts, 'sigma') && isstruct(opts.sigma)
    opts.tf_plot = opts.sigma;
end

%% Forward call to mess_tf_plot
opts.tf_plot.type = 'sigma';

[out, eqn, opts, oper] = mess_tf_plot(arg_one, opts, varargin{:});
