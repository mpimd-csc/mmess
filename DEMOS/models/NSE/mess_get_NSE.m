function [eqn, K0p, K0d] = mess_get_NSE(Re, level)
%% Load the linearized Navier Stokes example data and download it first,
% if it has not yet been downloaded.
%
% Call:
%  eqn = mess_get_NSE(Re, level)
%
% Inputs:
%  Re      desired Reynolds number (300, 400, 500)
%  level   desired refinement level
%            1 -> n =   3595 (velocity dim:  3142)
%            2 -> n =   9391 (velocity dim:  8268)
%            3 -> n =  22385 (velocity dim: 19770)
%            4 -> n =  50527 (velocity dim: 44744)
%            5 -> n = 110620 (velocity dim: 98054)
%
% Output:
%  eqn     equation structure for use with 'dae_2' usfs.
%  K0p     initial stabilizing feedback for the primal system
%  K0d     initial stabilizing feedback for the dual system

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%
opts = struct;
switch Re
    case {300, 400, 500}
    otherwise
        mess_err(opts, 'illegal_input', ...
                 'parameter ''Re'' must be 300, 400, or 500.');
end
base_file = sprintf('mat_nse_re_%i.mat', Re);
full_file = [fileparts(mfilename('fullpath')), filesep, base_file];
try
    load(full_file, 'mat');
catch
    mess_fprintf(opts, ['The file mat_nse_re_%i.mat, is used for the first time.', ...
                        'It is available as a separate download (270MB).\n\n'], Re);

    reply = input('Do you want to download it now? Y/N [Y]:', 's');
    if isempty(reply)
        reply = 'Y';
    end
    switch upper(reply)
        case 'Y'
            url = 'https://csc.mpi-magdeburg.mpg.de/';
            folder = 'mpcsc/software/mess/mmess/models/NSE/';
            mess_websave(full_file, [url, folder, base_file]);

        case 'N'
            mess_err(opts, 'check_data', ...
                     ['The download is required for NSE example. ', ...
                      'Consider switching to the simple Stokes model.']);
        otherwise
            mess_err(opts, 'illegal_input', 'Please answer Y or N.');
    end
    load(full_file, 'mat');
end

eqn.A_ = mat.mat_v.fullA{level};
eqn.E_ = mat.mat_v.E{level};
eqn.haveE = true;
eqn.B = mat.mat_v.B{level};
eqn.C = mat.mat_v.C{level};
eqn.manifold_dim = mat.mat_mg.nv(level);

K0p = mat.mat_v.Feed_0{level}';
K0d = mat.mat_v.Feed_1{level}';
if Re == 500
    K0d = K0d';
end
