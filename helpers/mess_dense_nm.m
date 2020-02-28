function X = mess_dense_nm(A,B,C,E,X0,S)
% naive Newton Kleinman iteration  for the ARE
%
%     A'*X*E+E'*X*A+C'*S*C-E'*X*B*B'*X*E = 0
%
% Inputs:
%
% A,B,C,E,S  Coefficients in the above equation
% X0         initial guess for the solution
%
% Outputs:
% X          Solution
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
%               2009-2020
%

% check for available Lyapunov solver
if exist('lyap', 'file')
    meth='lyap';
elseif exist('lyap_sgn_fac', 'file')
    % This should always be available since we ship it.
    % The following are only her n case someone deletes the sign solver
    meth='lyap_sgn_fac';
elseif exist('lyapchol', 'file')
    % In case lyap is not available it is rather unlikely that this will be
    % there. Still someone might have their own implementation...
    meth='lyapchol';
elseif exist('lyap2solve', 'file')
    % Again shipped with our code so actually unlikely to be unavailable.
    meth='lyap2solve';
else
    error('MESS:missing_solver',...
        'mess_dense_nm was unable to find Lyapunov solver');
end
tol = 1e-12;
maxiter = 50;

F = B*B';
if (nargin == 6) && not(isempty(S))
    if not(issymmetric(S))
        error('MESS:data', 'S must be symmetric');
    end
    G = C' * S * C;
    % S must be symmetric pos. semidef. for lyapchol or lyap_sgn_fac
    if strcmp(meth, 'lyapchol') || strcmp(meth, 'lyap_sgn_fac')
        [U,S_diag] = eig(S);
        if any(diag(S_diag)<0)
            meth='lyap2solve';
        end
    end
    switch meth
        case 'lyap'
            G = (G + G') / 2; % make sure it's symmetric for e.g. lyap
        case {'lyapchol', 'lyap_sgn_fac'}
            C = sqrt(S_diag) * U' * C; % C'*C = G
        case 'lyap2solve'
            if (nargin<4) || isempty(E)
                GE = C' * S * C;
            else
                CE = C/E;
                GE = CE'* S * CE;
            end
    end
else
    G = C' * C;
    switch meth
        case 'lyap'
            G = (G + G') / 2; % make sure it's symmetric for e.g. lyap
        case 'lyap2solve'
            if (nargin<4) || isempty(E)
                GE = C' * C;
            else
                CE = C/E;
                GE = CE' * CE;
            end
    end
end


res0 = norm(G);

if (nargin<4) || isempty(E)
    E = eye(size(A,1));
end
for i=1:maxiter
    if i>1 || ( (nargin==5) && not(isempty(X0)) )
        K = B'*X0*E;
    else
        K = zeros(size(B'));
        X0 = zeros(size(A));
    end
    switch meth
        case 'lyap'
            X = lyap(A'-K'*B',G+K'*K,[],E');
        case 'lyapchol'
            XC=lyapchol(A'-K'*B',[C',K'],E');
            X = XC'*XC;
        case 'lyap_sgn_fac'
            XC=lyap_sgn_fac(A - B*K,[C; K],E);
            X = XC'*XC;
        case 'lyap2solve'
            KE = K/E;
            X=lyap2solve(((A-B*K)/E)',GE+KE'*KE);
    end
    XE = X*E;
    res = norm(A'*XE+XE'*A-XE'*F*XE+G);
    rc = norm(X-X0)/norm(X);
    if (rc<tol) || (res<tol*res0)
        break
    else
        X0 = X;
    end
end
