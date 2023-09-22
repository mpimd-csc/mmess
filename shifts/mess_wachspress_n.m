function p = mess_wachspress_n(a, b, alpha, num_desired)
% function p=mess_wachspress_n(a,b,alpha,num_desired)
%
% calculates the optimal ADI shiftparameters (for equations where
% Matrix A is stable and symmetric) according to Jing-Rebecca Li
% and Jakob Whites "Low Rank Solution of Lyapunov equation" (which
% gives an overview of Wachspress`s method form e.g. "The ADI model Problem"
%
% p            is the array of shift parameters
%
% a            is assumed to be the absolute value of the smallest magnitude
%              eigenvalue of the Systemmatrix A
%
% b            is assumed to be the absolute value of the largest
%              magnitude eigenvalue
%
% alpha        is the arctan of the maximum of
%              abs(imag(lamba))/abs(real(lambda))
%              for all lambda in the spectrum of A
%
% num_desired  is the number of desired shift parameters that
%              should be calculated
%

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%
opts = struct;
if not(isnumeric(a)) || not(length(a) == 1)
    mess_err(opts, 'error_arguments', 'a has to be a numeric value');
end
if not(isnumeric(b)) || not(length(b) == 1)
    mess_err(opts, 'error_arguments', 'b has to be a numeric value');
end
if not(isnumeric(alpha)) || not(length(alpha) == 1)
    mess_err(opts, 'error_arguments', 'alpha has to be a numeric value');
end
if not(isnumeric(num_desired)) || not(length(num_desired) == 1)
    mess_err(opts, 'error_arguments', 'num_desired has to be a numeric value');
end

if alpha == 0
    num_Ritzrime = a / b;
else
    c2 = 2 / (1 + (a / b + b / a) / 2);
    m = 2 * cos(alpha) * cos(alpha) / c2 - 1;
    if m < 1
        mess_err(opts, 'complex_shifts', ...
                 ['Shift parameters would be complex, function not ' ...
                  'applicable, aborting!']);

        %
        % FIX ME: if m<1 parameters become complex! switch back to the
        % heuristics by Thilo or complex parameters suggested by
        % Wachspress.
        % If the reason are single outliers treat them separately (Wachspress
        % suggestion)
    end
    num_Ritzrime = 1 / (m + sqrt(m^2 - 1));
end
k = min(1 - eps, sqrt(1 - num_Ritzrime^2));
% this is a workaround for the case
% k=1 that works for Model reduction problems
% (not really nice but it works great for now).

% TODO: check the computation of k, kprime to avoid roundoff errors
% and probably replace the hack above.

[K, ~] = ellip(k, pi / 2);
J = num_desired;

p = ones(J, 1);
for i = 1:J
    p(i) = -sqrt(a * b / num_Ritzrime) * dn((i - 0.5) * K / J, k);
    % here we have the choice to take the
    % matlab function ellipj or our own
    % one dn. the later can be ported to
    % FORTRAN or C Code very easily
end
