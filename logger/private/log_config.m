function [out, format, figureformat, timeformat] = log_config
%% For every routine that does not get the 'opts' struct, these values will
%  be used as a default fallback to ensure the logger does not crash.

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

% these are the logging controls for every routine that does not get the
% 'opts' argument.
out          = 'console';
format       = 'md';
figureformat = 'png';
timeformat   = 'yyyy-MM-dd-HH-mm-ss';
