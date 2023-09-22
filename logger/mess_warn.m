function mess_warn(opts, reason, message, varargin)
%% warnings with additional output streams
% opts        the opts struct containing a logger field.
% reason      the warning code, restricted to the ones in mess_codes
% message     the warning message, describing the warning in finer
%             detail
% varargin    message can take sprintf-like arguments, these are the
%             arguments to sprintf(message,varargin)

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

% check for presence of opts.logger
opts = log_checkopts(opts);
% check for special case where warnings are turned on/off
% mess_warn(opts, 'OFF|ON', 'warning')

% check message dimensions
if not(size(message, 1) == 1)
    message = message';
    mess_warn(opts, 'warning_arguments', 'message should not be transposed');
end
disabled = warning('query');
disabledWarnings = struct2cell(disabled);
disabledWarnings = disabledWarnings(1, :);
if ismember(['MESS:' reason], disabledWarnings)
    return
end
if strcmpi(reason, 'ON') || strcmpi(reason, 'OFF')
    warning(reason, message);
    if not(strcmp(opts.logger.out, 'console'))
        warning('MESS:notimplemented', ...
                'warnings cannot be altered for file logging');
    end
    return
end

% varargin is not empty, replace placeholders in message
if not(isempty(varargin))
    message = sprintf(message, varargin{1:end});
end
codes = mess_log_codes('warning');

if not(sum(ismember(codes, reason) > 0))
    warning('MESS:error_arguments', ...
            ' %s is not a valid code', reason);
    return
end

file = [];
line = [];

db = dbstack(1);
if not(isempty(db))
    file = db(1).file;
    line = db(1).line;
end

MessageToPrint = log_formatter(opts, reason, message, 'warn');

% disable backtrace because that would show this file as a source, that's
% sort of wrong though.
warnstate = warning('off', 'backtrace');

switch opts.logger.out
    case 'file'

        fprintf(opts.logger.file, [MessageToPrint, '\n']);
    case 'both'
        fprintf(opts.logger.file, [MessageToPrint, '\n']);
        if not(isempty(file))
            warning(['MESS:' reason], ...
                    [message, '\n> in ', file, ' L', int2str(line)]);
            trace = interactiveStack(dbstack);
            for i = length(trace):-1:1
                fprintf(1, trace{i});
            end
        else % this is when mess_warn is called interactively
            warning(['MESS:' reason], message);
        end

    case 'console'
        if not(isempty(file))
            warning(['MESS:' reason], message);
            % print backtrace with opentoline refs
            trace = interactiveStack(dbstack);
            for i = length(trace):-1:1
                fprintf(1, trace{i});
            end
        else
            warning(['MESS:' reason], message);
        end

    otherwise
        error(MESS:illegal_log_location, ...
              'Requested unsupported log location');
end
% re-enable backtrace, so other stuff does not break
warning(warnstate);
