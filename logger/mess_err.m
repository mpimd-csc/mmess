function opts = mess_err(opts, reason, message, varargin)
%% error wrapper for additional output streams
% opts        the opts struct containing a logger field.
% reason      the error code, restricted to the ones in mess_codes
% message     the error message, describing the error in finer detail
% varargin    message can take sprintf-like arguments, these are the
%             arguments to sprintf(message,varargin)

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%
if not(size(message, 1) == 1)
    message = message';
    mess_warn(opts, 'warning_arguments', 'message should not be transposed');
end
% check for presence of opts.logger
opts = log_checkopts(opts);

% varargin is not empty, replace placeholders in message
if not(isempty(varargin))
    message = sprintf(message, varargin{1:end});
end

codes = mess_log_codes('error');
if not(sum(ismember(codes, reason) > 0))
    warning('MESS:error_arguments', 'error code is not valid');
    return
end

db = dbstack();

MessageToPrint = log_formatter(opts, reason, message, 'err');

switch opts.logger.out
    case 'console'
        % nothing to do. the actual error is handled below.

    case 'file'
        fprintf(opts.logger.file, [MessageToPrint, '\n']);
        mess_log_finalize(opts);

    case 'both'
        fprintf(opts.logger.file, [MessageToPrint, '\n']);
        opts = mess_log_finalize(opts);

    otherwise
        error(MESS:illegal_log_location, ...
              'Requested unsupported log location');

end

% the program should halt under any circumstances when an error is thrown,
% otherwise this would not qualify as an error

% make MATLAB throw the error three levels above this file,
% this way not every error is thrown in mess_err
errorStruct.message = char(message);
errorStruct.identifier = char(strcat('MESS:', reason));
errorStruct.stack = db(2:end);
% throw the error
error(errorStruct);

end
