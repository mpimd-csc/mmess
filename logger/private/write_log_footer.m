function write_log_footer(opts)
%%  writes the matching file header for the specified format
%   format      specifies which header is to be returned
%
%   header      the first lines to print into the file for <format>-logging

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%
file = ['log.' lower(opts.logger.format) '.templ'];

% find the templates
messpath = fileparts(which('mess_path.m'));
template = [messpath filesep 'logger' filesep 'resources' filesep file];
founddiv = false;
fid = fopen(template);
tline = fgetl(fid);
if strcmp(tline, '======HEADER ABOVE, FOOTER BELOW, DO NOT MODIFY THIS LINE======')
    founddiv = true;
end
while ischar(tline)
    if founddiv
        fprintf(opts.logger.file, '%s\n', tline);
    end

    if strcmp(tline, '======HEADER ABOVE, FOOTER BELOW, DO NOT MODIFY THIS LINE======') && not(founddiv)
        founddiv = true;
    end
    tline = fgetl(fid);
end
fclose(fid);
