function mess_log_plot(opts, figure, figurename)
%% logging for figures into files as well as embedding into text files
% opts        the opts struct containing a logger field.
% figure      figure handle to the figure
% figurename  string naming the figure

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

if nargin < 3
    % octave compatibility
    if not(exist('OCTAVE_VERSION', 'builtin'))
        figurename = strcat('figure', string(figure.Number), '.png');
    else
        figurename = strcat('figure', num2str(figure), '.png');
    end
end
supported_formats = {'png', 'eps', 'svg'};

opts = log_checkopts(opts);

if isfield(opts.logger, 'figureformat') && ...
        ismember(opts.logger.figureformat, supported_formats)

    type = opts.logger.figureformat;
else
    type = 'png';
end

if not(strcmp(opts.logger.out, 'console'))
    full_figurename = [opts.logger.messlogdir, filesep, figurename];
    saveas(figure, full_figurename, type);
    saveas(figure, full_figurename, 'png');

    switch opts.logger.format
        case 'md'
            fprintf(opts.logger.file, ...
                    '![%s](%s) \n', ...
                    figurename, ...
                    [full_figurename, '.', type]);
        case 'tex'
            fprintf(opts.logger.file, ...
                    '\\includegraphics{%s} \\\\ \n', ...
                    [figurename, '.', type]);
        case 'html'
            fprintf(opts.logger.file, ...
                    '<img src="%s"><br> ', ...
                    [full_figurename, '.', type]);
        otherwise
            % error treating
    end
end
