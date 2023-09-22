function opts = mess_log_finalize(opts)
%% finalizes the logging, renders pdf in some cases
% opts        the opts struct containing a logger field.

% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)

opts = log_checkopts(opts);

if strcmp(opts.logger.out, 'console')
    return
end
write_log_footer(opts);

fclose(opts.logger.file);
opts.logger = rmfield(opts.logger, 'file');
