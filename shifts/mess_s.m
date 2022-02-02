function [max_r,ind] = mess_s(p,set)
%
% Computation of the maximal magnitude of the rational ADI function over
% a discrete subset of the left complex half plane.
%
%   Calling sequence:
%
%     [max_r,ind] = mess_s(p,set)
%
%   Input:
%
%     p        vector of ADI parameters;
%     set      vector representing the discrete set.
%
%   Output:
%
%     max_r    maximal magnitude of the rational ADI function over set;
%     ind      index - maximum is attained for set(ind).
%

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%


%   Exact copy from
%
%   LYAPACK 1.0 (Thilo Penzl, Jan 1999)
%
if not(isnumeric(p))
    error('MESS:error_arguments','p has to be a vector of numeric type');
end
if not(isnumeric(set))
    error('MESS:error_arguments','set has to be a vector of numeric type');
end
max_r = -1;
ind = 0;

for i = 1:length(set)

  x = set(i);

  rr = 1;
  for j = 1:length(p)

    rr = rr*abs(p(j)-x)/abs(p(j)+x);

  end

  if rr > max_r

    max_r = rr;
    ind = i;

  end

end



