function opts = log_checkopts(opts)
%% The opts will be checked for completeness, if opts is (partially)
%  missing, it will be filled with the default values from "log_config"
%
%  opts             potentially incomplete options struct, concerning the
%                   opts.logger part
%
%
%  opts             logging-wise complete options struct.

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

% default options can be configured in log_config.m
[out, format, figureformat, timeformat] = log_config();

% checking opts, and defaulting to console logging
if not(isfield(opts, 'logger'))
    opts.logger = struct();
end
if not(isfield(opts.logger, 'out'))
    opts.logger.out = out;
end
% in case someone forgot the format, default to md
if not(isfield(opts.logger, 'format'))
    opts.logger.format = format;
end
if not(isfield(opts.logger, 'figureformat'))
    opts.logger.figureformat = figureformat;
end
if not(isfield(opts.logger, 'basename'))
    if exist('OCTAVE_VERSION', 'builtin')
        opts.logger.basename = datestr(now(), timeformat);
    else
        opts.logger.basename = char(datetime('now', 'Format', timeformat));
    end
end
