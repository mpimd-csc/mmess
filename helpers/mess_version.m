function [ver, opts] = mess_version(opts, quiet)
% MESS_VERSION prints a short version message and returns the version
% number as a string
%
% Input
%
%  opts  options structure containing the logger data for correct
%        printing of the message
%        (optional, defaults to empty struct, i.e. printing only to
%         the console)
%  quiet switch printing of message of when set to 'quiet'
%        (optional, default is printing enabled)
%
% Output
%
%  ver   the numeric version as a string
%  opts  the options structure
%

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

ver = '3.0';

if nargin < 2 || not(strcmp(quiet, 'quiet'))
    if nargin < 1 || isempty(opts)
        opts = mess_log_initialize(struct());
    else
        if not(isfield(opts, 'log'))
            mess_log_initialize(opts);
        end
    end

    mess_fprintf(opts, '\n');
    mess_fprintf(opts, 'This is M-M.E.S.S. version %s\n\n', ver);
end

end
