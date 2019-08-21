function p = mess_mnmx(rw,num_desired)
%
%  Suboptimal solution of the ADI minimax problem. The delivered parameter
%  set is closed under complex conjugation. 
%
%  Calling sequence:
%
%    p = mess_mnmx(rw,num_desired)
%
%  Input:
%
%    rw            a vector containing numbers in the open left
%                  half plane, which approximate the spectrum of
%                  the corresponding matrix, e.g., a set of Ritz
%                  values. The set must be closed w.r.t. complex
%                  conjugation; 
%    num_desired   desired number of shift parameters 
%                  (length(rw) >= num_desired)
%                  (The algorithm delivers either num_desired or
%                  num_desired+1 parameters, to ensure closedness
%                  unde complex conjugation!).
%
%  Output:
%
%    p             a num_desired- or num_desired+1-vector of
%                  suboptimal ADI parameters;
%

%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, see <http://www.gnu.org/licenses/>.
%
% Copyright (C) Jens Saak, Martin Koehler, Peter Benner and others 
%               2009-2019
%

%  Exact copy from
%   
%  LYAPACK 1.0 (Thilo Penzl, October 1999)

% Input data not completely checked!
if not(isnumeric(num_desired)) || (length(num_desired) ~= 1)
    error('MESS:error_arguments','num_desired has to be of numeric type');
end
if not(isnumeric(rw))
    error('MESS:error_arguments','rw has to be a vector of numeric type');
end
if length(rw)<num_desired
  error('length(rw) must be at least num_desired.');
end

max_rr = +Inf;                       % Choose initial parameter (pair)

for i = 1:length(rw)
  max_r = mess_s(rw(i),rw); 
  if max_r < max_rr
    p0 = rw(i);
    max_rr = max_r;
  end
end  

if imag(p0)
  p = [ p0; conj(p0) ];
else
  p = p0;                            
end  

[~,i] = mess_s(p,rw);         % Choose further parameters.

while size(p,1) < num_desired 
   
  p0 = rw(i);
  if imag(p0)
    p = [ p; p0; conj(p0) ]; %#ok<AGROW>
  else
    p = [ p; p0]; %#ok<AGROW>
  end
  
  [~,i] = mess_s(p,rw);
    
end

