function mess_fprintf(opts, message, varargin)
%% prints the message to the output(s) specified in opts
% opts        the opts struct containing a logger field.
% message     the message, can contain format specifying placeholders
% varargin    sprintf-like arguments to 'message'

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

% check for presence of opts.logger
opts = log_checkopts(opts);
if not(size(message, 1) == 1)
    message = message';
    mess_warn(opts, 'warning_arguments', 'message should not be transposed');
end
% varargin is not empty, replace placeholders in message
if not(isempty(varargin))
    if strcmp(opts.logger.format, 'md')
        message =  strrep(message, '\n', '  \n');
    end
    message = sprintf(message, varargin{1:end});
end

MessageToPrint = log_formatter(opts, [], message, 'log');

switch opts.logger.out
    case 'console'
        fprintf(message);

    case 'file'
        if strcmp(opts.logger.format, 'md')
            fprintf(opts.logger.file, MessageToPrint);
        else
            fprintf(opts.logger.file, [MessageToPrint, '\n']);
        end
    case 'both'
        if strcmp(opts.logger.format, 'md')
            fprintf(opts.logger.file, MessageToPrint);
        else
            fprintf(opts.logger.file, [MessageToPrint, '\n']);
        end
        fprintf(message);

    otherwise
        error(MESS:illegal_log_location, ...
              'Requested unsupported log location');

end
