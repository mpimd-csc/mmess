function data = mess_load_rail(k)
%% function data = mess_load_rail(k) reads in a set of normalized
% finite element data matrices for the rail example
% s.a. Oberwolfach Model Reduction Benchmark Collection hosted at MORWiki
% (https://modelreduction.org)
%
% DESCRIPTION:
%   This function extracts the mat-files to get the linear steel rail
%   model described by
%
%       E*x' = A*x + B*u                                               (1a)
%          y = C*x,                                                    (1b)
% 
%   or the bilinear steel rail model described by
%
%       E*x' = A*x + B*u + N{1}*x*u_1 + ... + N{6}*x*u_6,              (1a)
%          y = C*x,                                                    (1b)
%
%   for different refinement sizes. This function uses the FENICS
%   reimplementation of the Oberwolfach Collection steel profile
%   and just like the original model features 7 inputs and 6
%   outputs. The data is intended for post-processing 
%
% Input:
%   k       number of instance
%           k = 0 -> n =       109
%           k = 1 -> n =       371
%           k = 2 -> n =     1,357
%           k = 3 -> n =     5,177
%           k = 4 -> n =    20,209
%           k = 5 -> n =    79,841
%           k = 6 -> n =   317,377
%           k = 7 -> n = 1,265,537
%           k = 8 -> n = 5,054,209
%
% Output:
%  data  structure containing the separate normalized matrices used
%        to generate the model variants in the post processing.
%     
%        members are:
%         M  the mass matrix
%         S  the stiffness matrix (aka discretized Laplace)
%         M_Gamma_k for k = 0, ... , 6 the boundary mass matrices
%         B_k       for k = 0, ... , 6 the discretized spatial
%                   input distributions on the corresponding boundaries.
%
% Data generated by:
%   https://gitlab.mpi-magdeburg.mpg.de/models/fenicsrail
%
% REFERENCES:
%  [1] J. Saak, Efficient numerical solution of large scale
%      algebraic matrix equations in PDE control and model order
%      reduction, Dissertation, Technische Universität Chemnitz,
%      Chemnitz, Germany (Jul. 2009).
%      URL http://nbn-resolving.de/urn:nbn:de:bsz:ch1-200901642
%  [2] J. Saak, Effiziente numerische Lösung eines
%      Optimalsteuerungsproblems für die Abkühlung von Stahlprofilen,
%      Diplomarbeit, Fachbereich 3/Mathematik und Informatik,
%      Universität Bremen, D-28334 Bremen (Sep. 2003).
%  [3] P. Benner, J. Saak, A semi-discretized heat transfer model
%      for optimal cooling of steel profiles, in: P. Benner,
%      V. Mehrmann, D. Sorensen (Eds.), Dimension Reduction of
%      Large-Scale Systems, Vol. 45 of Lect. Notes Comput. Sci. Eng.,
%      Springer-Verlag, Berlin/Heidelberg, Germany, 2005,
%      pp. 353–356. https://doi.org/10.1007/3-540-27909-1_19.

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%% set path
switch k
    case 0
        example = 'ODE_unit_matrices_109';
    case 1
        example = 'ODE_unit_matrices_371';
    case 2
        example = 'ODE_unit_matrices_1357';
    case 3
        example = 'ODE_unit_matrices_5177';
    case 4
        example = 'ODE_unit_matrices_20209';
    case 5
        example = 'ODE_unit_matrices_79841';
    case 6
        example = 'ODE_unit_matrices_317377';
    case 7
        example = 'ODE_unit_matrices_1265537';
    case 8
        example = 'ODE_unit_matrices_5054209';
    otherwise
        error('MESS:error_arguments', ...
              'k must be a non-negative integer smaller than 9.');
end

%% check path
path = fileparts(mfilename('fullpath'));
base_file = [example '.mat'];
full_file = [path filesep '..' filesep base_file];

%% read matrices
try
    data = load(full_file);
catch
    % file not (yet?) available (we are only shipping the small
    % resolutions to keep the download size limited)
    fprintf(['The file %s, is used for the first time. ' ...
        'It is available as a separate download.\n\n'], base_file);

    % Let's warn for the really large dowload sizes
    switch k
        case 7
            fprintf('Attention: This file is 329M large!\n');
        case 8
            fprintf('Attention: This file is 1.3G large!\n');
    end
    reply = input('Do you want to download it now? Y/N [Y]:','s');

    if isempty(reply)
        reply = 'Y';
    end

    switch upper(reply)
        case 'Y'
            % Download requested. This is for the location:
            url = 'https://csc.mpi-magdeburg.mpg.de/';
            folder = 'mpcsc/software/mess/mmess/models/Rail/';
            
            % If possible use the modern websave funciton otherwise
            % (especially in Octave) fall back to the classic
            % urlwrite to fetch the matrix data from the above location
            if exist('websave','file')
                websave(full_file, [url, folder, base_file]);
            else
                urlwrite([url, folder, base_file], full_file); %#ok<URLWR>
            end
            
        case 'N'
            % we can not continue without the file.         
            error('The download is required.');
            
        otherwise
            % illegal user input. 
            error('Please answer Y or N.');
    end
    
    % OK, we should have the file now, let finally retry to load it.
    data = load(full_file);
    
end

end
