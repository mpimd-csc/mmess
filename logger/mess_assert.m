function opts = mess_assert(opts, condition, reason, message, varargin)
%% assert wrapper for additional output streams
% opts        the opts struct containing a logger field.
% condition   condition to assert
% reason      the error code to throw if the assertion fails,
%             restricted to the ones in mess_codes
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

if not(condition)
    mess_err(opts, reason, message, varargin);
end

end
