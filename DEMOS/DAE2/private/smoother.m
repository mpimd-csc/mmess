function alpha = smoother(alpha, i, ntau)
% Smoother for implicit euler

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

if i < (0.2 * ntau)

    alpha = sin((10.0 * pi * (i - 0.1 * ntau) / ntau) - 0.5 * pi) * alpha;
end
end
