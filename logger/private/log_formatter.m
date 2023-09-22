function formatted_Message = log_formatter(opts, reason, message, level)
%% formats the reason and message according to the opts and the level
% opts              the options struct
% reason            the reason the exception was thrown
% message           further explanation on the error
% level             corresponds to the calling function;
%                   'log'                fprintf
%                   'warn'               warn
%                   'err'                error
%
% formatted_Message combined reason & message according to opts & level

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

stack = dbstack;
% stack(1) is this,
% stack(2) is mess_err/mess_warn
% stack(3) is whatever calls mess_err/mess_warn
if length(stack) > 2
    line = stack(3).line;
    file = stack(3).file;
else
    line = 0;
    file = 'command-line';
end

delim1 = '';
delim2 = '';
IDstring = 'MESS:';
mark = '>';

switch opts.logger.format
    case 'tex'
        newline  = '\\newline{}';
        mark = '$>$';
        file = strrep(file, '_', '\\_');
        switch level
            case 'warn'
                delim1 = '\\textcolor{orange}{';
                delim2 = '}';
                keyword = '-Warning: ';
            case 'err'
                delim1 = '\\textcolor{red}{';
                delim2 = '}';
                keyword = '-Error: ';
        end

        % the reason mostly contains stuff like <illegal_input>
        % to escape this properly to latex this needs to become
        % <illegal\\_input> because fprintf will choke otherwise
        if not(isempty(reason))
            reason = cleantex(reason);
        end

        message = cleantex(message);
    case 'md'
        newline = '';
        switch level
            case 'warn'
                delim1 = '*';
                delim2 = '*';
                keyword = '-Warning: ';
            case 'err'
                delim1 = '***';
                delim2 = '***';
                keyword = '-Error: ';
        end
        message = strrep(message, '\n', newline);
    case 'html'
        newline = '<br/>';
        switch level
            case 'warn'
                delim1 = '<p style="color:rgb(255,146,0);">';
                delim2 = '</p>';
                keyword = '-Warning: ';
            case 'err'
                delim1 = '<p style="color:rgb(255,0,0);">';
                delim2 = '</p>';
                keyword = '-Error: ';
        end
        message = strrep(message, '\n', newline);
end
if not(isempty(reason))
    reason = [reason, keyword];
    if strcmp(opts.logger.format, 'md')
        formatted_Message = ...
          [delim1, IDstring, reason, message, newline, mark, ' in ', ...
           file, ' L', int2str(line), delim2, newline];
    else
        formatted_Message = ...
          [delim1, ' ', IDstring, reason, message, newline, mark, ' in ', ...
           file, ' L', int2str(line), delim2, newline];
    end
else % 'log' level case
    if strcmp(opts.logger.format, 'md')
        formatted_Message = message;
    else
        formatted_Message = [message, newline];
    end
end
