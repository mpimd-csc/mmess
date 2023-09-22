function logging(setup, is_test)
% Demo for logging functions of M-M.E.S.S.

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

% set the options for the logger in the opts.logger struct
% figureformat sets the format of the graphs/figures that will be used
% (eps, png, svg)
opts.logger.figureformat = 'png';

% format sets the format of the document that will be generated
% (html, md, tex)
opts.logger.format       = 'md';

% out sets the output:
%   'file'     redirect output to file,
%   'console'  print to the matlab console as usual
%   'both'     print as usual but duplicate the output into the file
opts.logger.out          = 'file';

% fetch the above from 'setup' input argument in case this is executed
% as a unit test
if nargin > 0
    opts.logger.format       = setup{1};
    opts.logger.out          = setup{2};
    opts.logger.figureformat = setup{3};
    opts = mess_log_initialize(opts, ...
                               [setup{3} '-' setup{1} '-' setup{2}]);
end

if nargin < 2
    is_test = false;
end

if not(is_test)
    % initialize the logging mechanisms by giving it a name
    % this needs to be written to the opts struct and called after the
    % setting of the logger options
    opts = mess_log_initialize(opts, 'demo_log');
end

%% mess_fprintf works like the regular fprintf, but will print to the output
% set in opts.logger.out
mess_fprintf(opts, ...
             '============ %s =============\n', ...
             opts.logger.basename);

mess_fprintf(opts, ...
             '%s %s\n', ...
             opts.logger.format, opts.logger.out);

% print some numbers
for j = 1:10
    mess_fprintf(opts, 'test: %2d\n', j);
end

%% mess_warn is used to issue warnings, the warning codes are restricted to
% a set that can be found in mess_log_codes.m
mess_warn(opts, ...
          'illegal_input', ...
          'this warning was triggered by incorrect inputs');

%% logging plots
% plot a quadratic
x = linspace(-5, 5);
y = x.^2;
f = figure();
plot(x, y);
xlabel('x axis');
ylabel('y axis');
title('demo figure');

% mess_log_plot will print the given figure in the format set in
% opts.logger.figureformat. the file will be called NAME.{figureformat}
% if NAME is not given, the figure's number is used.
NAME = 'demo';
mess_log_plot(opts, f, NAME);
close(f);

%% log a matrix
m = sprand(10, 10, 0.1);
% mess_log_matrix saves the given variable to a .mat file
mess_log_matrix(opts, m);

if is_test
    try
        %% mess_err throws an error, logs to the output set in opts.logger.out
        % and properly closes the file if necessary.
        mess_err(opts, ...
                 'missing_feature', ...
                 'this error was triggered by a missing feature');
    catch

    end

else
    %% mess_log_finalize needs to be called at the end of every logged
    % computation, it finishes up the logging process.
    % note that this is also automatically called in mess_err, should the
    % computation end due to an error
    mess_log_finalize(opts);
end
