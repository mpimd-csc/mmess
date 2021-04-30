%  Driver script for the results presented in:
%
%  J. Saak, M. Voigt, Model reduction of constrained mechanical systems
%  in M-M.E.S.S., IFAC-PapersOnLine 9th Vienna International Conference
%  on Mathematical Modelling MATHMOD 2018, Vienna, Austria, 21–23
%  February 2018 51 (2) (2018) 661–666.
%  https://doi.org/10.1016/j.ifacol.2018.03.112

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2021 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

model={{'Stykel_large', 1e-4, 200, 500}, ...
       {'Truhar_Veselic', 1e-4, 250, 300}, ...
       {'TV2', 1e-3, 300, 300},...
       {'TV2', 1e-3, 300, 500}...
      };
l = length(model);
for i=1:l
    BT_DAE3_SO(model{i}{1},model{i}{2},model{i}{3},model{i}{4})
    if (i<l)
        fprintf('\nPress any key to continue\n');
        pause; 
    end
end
