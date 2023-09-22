function opts = mess_log_initialize(opts, name)
%% initializes the logging mechanisms
%
% opts        the opts struct containing a logger field.
%
% name        the working title of this computation, used as the basename
%             for the log-file.
%             (optional, defaults to name of the calling function)
%
% returns:    the updated opts struct

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)

if not(isfield(opts, 'logger')) || not(isfield(opts.logger, 'out'))
    opts = log_checkopts(opts);
end

if not(strcmp(opts.logger.out, 'console'))
    if exist('name', 'var')
        opts.logger.basename = name;
    else
        ST =  dbstack;
        if length(ST) > 1
            opts.logger.basename = ST(2).name;
            name =  ST(2).name;
        else
            name = 'interactive';
        end
    end

    % create the ./mess_log folder, if not present
    logpath = pwd;
    [~, ~, ~, timeformat] = log_config;
    if exist('OCTAVE_VERSION', 'builtin')
        timestamp = datestr(now(), timeformat);
    else
        timestamp = char(datetime('now', 'Format', timeformat));
    end
    logdir = ['mess_log-', name, '-', timestamp];

    [SUCCESS, MESSAGE, ~] = mkdir(logpath, logdir);
    if not(SUCCESS)
        error('MESS:logger_init', ...
              'Could not create M.E.S.S. log directory .\n Reason given: %s', ...
              MESSAGE);
    end

    % opts must be complete before this is called
    opts = log_checkopts(opts);

    opts.logger.messlogdir = [logpath, filesep, logdir];
    filename = [opts.logger.basename, '.', opts.logger.format];

    [opts.logger.file, errmsg] = ...
        fopen([opts.logger.messlogdir, filesep, filename], 'w+');

    if opts.logger.file < 0
        error('MESS:logger_init', ...
              'Could not create M.E.S.S. log file.\n Reason given: %s', ...
              errmsg);
    end

    write_log_header(opts);

end
