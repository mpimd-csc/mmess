function [out, eqn, opts, oper] = ...
    mess_Frobenius_TF_error_plot(arg_one, opts, varargin)
% Computation of simple Frobenius-norm-magnitude-plots for descriptor
% systems and comparison to reduced order models.
%
% Backward compatibility wrapper for mess_tf_plot with
%    opts.tf_plot.type == 'frobenius'
%
% arguments are as in mess_tf_plot. All settings formerly passed as
% arguments, now need to be put to opts.tf_plot.
%

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%% Forward call to mess_tf_plot
opts.tf_plot.type = 'frobenius';

[out, eqn, opts, oper] = mess_tf_plot(arg_one, opts, varargin{:});
