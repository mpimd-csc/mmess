function eqn = mess_get_BIPS(model, alpha, opts)
%% Load the power systems example data and download it first,
% if it has not yet been downloaded.
% Examples taken from https://sites.google.com/site/rommes/software
%
% Call:
%  eqn = mess_get_BIPS(model, alpha)
%
% Inputs:
%  model      BIPS model selection (integer, default: 7)
%                name           size  inputs outputs
%              1 BIPS/97        13251   1      1
%              2 BIPS/1997      13250   1      1
%              3 BIPS/97        13309   8      8
%              4 BIPS/97        13251  28     28
%              5 BIPS/97        13250  46     46
%              6 Juba5723       40337   2      1
%              7 bips98_606      7135   4      4
%              8 bips98_1142     9735   4      4
%              9 bips98_1450    11305   4      4
%             10 bips07_1693    13275   4      4
%             11 bips07_1998    15066   4      4
%             12 bips07_2476    16861   4      4
%             13 bips07_3078    21128   4      4
%  alpha     positive and real alpha in the alpha shift strategy
%            suggested by Rommes and coauthors in [1].
%            (optional, default: 0.05)
%  opts      M.E.S.S. options for the logger
%
%
% Output:
%  eqn     equation structure for use with 'dae_1' usfs.
%
% References:
% [1] F. Freitas, J. Rommes, N. Martins, Gramian-based reduction method
%     applied to large sparse power system descriptor models, IEEE Transactions
%     on Power Systems 23 (3) (2008) 1258–1270. doi:10.1109/TPWRS.2008.926693
%

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

if nargin < 1
    model = 7;
end

if nargin < 2
    alpha = 0.05;
end

if nargin < 3
    opts = struct();
end

filenames = { ...
             'siso_ww_vref_6405.mat', ...
             'siso_xingo_afonso_itaipu.mat', ...
             'mimo8x8_system.mat', ...
             'mimo28x28_system.mat', ...
             'mimo46x46_system.mat', ...
             'juba40k.mat', ...
             'bips98_606.mat', ...
             'bips98_1142.mat', ...
             'bips98_1450.mat', ...
             'bips07_1693.mat', ...
             'bips07_1998.mat', ...
             'bips07_2476.mat', ...
             'bips07_3078.mat'
            };

base_file = filenames{model};
full_file = [fileparts(mfilename('fullpath')), filesep, base_file];
try
    Bips = load(full_file);
catch
    mess_fprintf(opts, ['The file %s, is used for the first time.', ...
                        'It is available as a separate download.\n\n'], base_file);

    reply = input('Do you want to download it now? Y/N [Y]:', 's');

    if isempty(reply)
        reply = 'Y';
    end

    switch upper(reply)
        case 'Y'
            url = 'https://csc.mpi-magdeburg.mpg.de/';
            folder = 'mpcsc/software/mess/mmess/models/BIPS/';
            mess_websave(full_file, [url, folder, base_file]);

        case 'N'
            mess_err(opts, 'check_data', 'The download is required.');
        otherwise
            mess_err(opts, 'illegal_input', 'Please answer Y or N.');
    end
    Bips = load(full_file);
end
p = find(diag(Bips.E));
np = find(diag(Bips.E) == 0);
pp = [p; np];
if alpha == 0
    eqn.A_ = Bips.A(pp, pp);
else
    eqn.A_ = Bips.A(pp, pp) - alpha * Bips.E(pp, pp);
end
eqn.E_ = Bips.E(pp, pp);
switch model
    case {7, 8, 9, 10, 11, 12, 13}
        eqn.B = Bips.b(pp, :);
        eqn.C = Bips.c(:, pp);
    case {1, 3, 4}
        eqn.B = Bips.b(pp, :);
        eqn.C = Bips.c(pp, :)';
    case {2, 5}
        eqn.B = Bips.B(pp, :);
        eqn.C = Bips.C(pp, :)';
    case 6
        eqn.B = Bips.B(pp, :);
        eqn.C = Bips.C(:, pp);
        mess_fprintf(opts, 'Attention: This model is not stable!');
end
eqn.manifold_dim = length(p);
eqn.haveE = true;
