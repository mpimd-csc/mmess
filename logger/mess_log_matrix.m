function mess_log_matrix(opts, varargin)
%% logging matrices into files
% opts        the opts struct containing a logger field.
% varargin    data that should be logged, can be multiple variables

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

opts = log_checkopts(opts);
if nargin
    % A matrix is logged.
    for i = 1:length(varargin)
        filename = [opts.logger.messlogdir, filesep, ...
                    inputname(i + 1), '.mat'];
        assign(inputname(i + 1), varargin(i));
        save(filename, inputname(i + 1));
        if isfield(opts, 'logger') && isfield(opts.logger, 'out') && ...
                not(strcmp(opts.logger.out, 'console'))
            mess_fprintf(opts, ...
                         '%s was logged to %s \n', ...
                         inputname(i + 1), ...
                         filename);
        end
    end
end
end

function assign(VarName, VarValue)
assignin('caller', VarName, VarValue);
end
