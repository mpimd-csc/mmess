function C = mul_ApE_dae_1(eqn, opts, opA, p, opE, B, opB)%#ok<INUSL>

%% function mul_A perfoms operation C = (opA(A_)+pc*opE(E_))*opB(B)
%
% Input:
%   eqn     structure contains A_
%
%   opts    struct contains parameters for the algorithm
%
%   opA     character specifies the form of opA(A_)
%           opA = 'N' performs A_*opB(B)
%           opA = 'T' performs A_'*opB(B)
%
%   B       m-x-p matrix
%
%   opB     character specifies the form of opB(B)
%           opB = 'N' performs opA(A_)*B
%           opB = 'T' performs opA(A_)*B'
%
% Output:
% C = (opA(A_)+pc*opE(E_))*opB(B)
%
%   uses size_dae_1

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

%% check input Paramters
if (not(ischar(opA)) || not(ischar(opB)) || not(ischar(opE)))
    error('MESS:error_arguments', 'opA, opB or opE is not a char');
end

opA = upper(opA); opB = upper(opB); opE = upper(opE);

if(not((opA == 'N' || opA == 'T')))
    error('MESS:error_arguments','opA is not ''N'' or ''T''');
end

if(not((opB == 'N' || opB == 'T')))
    error('MESS:error_arguments','opB is not ''N'' or ''T''');
end

if(not((opE == 'N' || opE == 'T')))
    error('MESS:error_arguments','opE is not ''N'' or ''T''');
end

if(not(isnumeric(p)))
    error('MESS:error_arguments','p is not numeric');
end

if not(isfield(eqn, 'haveE')), eqn.haveE = 0; end

if (not(isnumeric(B))) || (not(ismatrix(B)))
    error('MESS:error_arguments','B has to ba a matrix');
end

%% check data in eqn structure
if(eqn.haveE ==1)
    if(not(isfield(eqn,'E_')) || not(isnumeric(eqn.E_))...
            || not(isfield(eqn,'A_'))) || not(isnumeric(eqn.A_))
        error('MESS:error_arguments','field eqn.E_ or eqn.A_ is not defined or corrupted');
    end
else
    if(not(isfield(eqn,'A_'))) || not(isnumeric(eqn.A_))
        error('MESS:error_arguments','field eqn.A_ is not defined');
    end
end
if not(isfield(eqn, 'st'))    || not(isnumeric(eqn.st))
    error('MESS:st',...
        'Missing or Corrupted st field detected in equation structure.');
end
n = size(eqn.A_,1);
st = eqn.st;
%% perfom multiplication
switch opA
    
    case 'N'
        switch opE
            case 'N'
                switch opB
                    
                    %implement operation A_*B
                    case 'N'
                        if(st > size(B,1))
                            error('MESS:error_arguments', 'number of cols of A_ differs with rows of B');
                        end
                        C = (eqn.A_(1 : st, 1 : st) + p * eqn.E_(1 : st, 1 : st)) * B(1 : st, : ) ...
                            - eqn.A_(1 : st, st + 1 : n) ...
                            * (eqn.A_(st + 1 : n, st + 1 : n) \ (eqn.A_(st + 1 : n, 1 : st) ...
                            * B(1 : st, : )));
                        
                        %implement operation A_*B'
                    case 'T'
                        if(st > size(B, 2))
                            error('MESS:error_arguments', 'number of cols of A_ differs with cols of B');
                        end
                        C = (eqn.A_(1 : st, 1 : st) + p * eqn.E_(1 : st, 1 : st)) * B( : , 1 : st)' ...
                            - eqn.A_(1 : st, st + 1 : n) ...
                            * (eqn.A_(st + 1 : n, st + 1 : n) \ (eqn.A_(st + 1 : n, 1 : st) ...
                            * B( : , 1 : st)'));
                end
            case 'T'
                switch opB
                    
                    %implement operation A_*B
                    case 'N'
                        if(st > size(B,1))
                            error('MESS:error_arguments', 'number of cols of A_ differs with rows of B');
                        end
                        C = (eqn.A_(1 : st, 1 : st) + p * eqn.E_(1 : st, 1 : st)') * B(1 : st, : ) ...
                            - eqn.A_(1 : st, st + 1 : n) * ...
                            (eqn.A_(st + 1 : n, st + 1 : n) \ (eqn.A_(st + 1 : n, 1 : st) * ...
                            B(1 : st, : )));
                        
                        %implement operation A_*B'
                    case 'T'
                        if(st > size(B, 2))
                            error('MESS:error_arguments', 'number of cols of A_ differs with cols of B');
                        end
                        C = (eqn.A_(1 : st, 1 : st) + p * eqn.E_(1 : st, 1 : st)') * B( : , 1 : st)'...
                            - eqn.A_(1 : st, st + 1 : n) ...
                            * (eqn.A_(st + 1 : n, st + 1 : n) \ (eqn.A_(st + 1 : n, 1 : st) ...
                            * B( : , 1 : st)'));
                end
        end
        
    case 'T'
        switch opE
            case 'N'
                switch opB
                    
                    %implement operation A_'*B
                    case 'N'
                        if(st > size(B, 1))
                            error('MESS:error_arguments', 'number of rows of A_ differs with rows of B');
                        end
                        C = (eqn.A_(1 : st, 1 : st)' + p * eqn.E_(1 : st, 1 : st)) * B(1 : st, : ) ...
                            - eqn.A_(st + 1 : n, 1 : st)' ...
                            * (eqn.A_(st + 1 : n, st + 1 : n)' \ (eqn.A_(1 : st, st + 1 : n)' ...
                            * B(1 : st, : )));
                        
                        %implement operatio A_'*B'
                    case 'T'
                        if(st > size(B, 2))
                            error('MESS:error_arguments', 'number of rows of A_ differs with cols of B');
                        end
                        C = (eqn.A_(1 : st, 1 : st)' + p * eqn.E_(1 : st, 1 : st)) * B( : , 1 : st)'...
                            - eqn.A_(st + 1 : n, 1 : st)' ...
                            * (eqn.A_(st + 1 : n, st + 1 : n)' \ (eqn.A_(1 : st, st + 1 : n)' ...
                            * B( : , 1 : st)'));
                end
            case 'T'
                switch opB
                    
                    %implement operation A_'*B
                    case 'N'
                        if(st > size(B, 1))
                            error('MESS:error_arguments', 'number of rows of A_ differs with rows of B');
                        end
                        C = (eqn.A_(1 : st, 1 : st)' + p * eqn.E_(1 : st, 1 : st)') * B(1 : st, : ) ...
                            - eqn.A_(st + 1 : n, 1 : st)' ...
                            * (eqn.A_(st + 1 : n, st + 1 : n)' \ (eqn.A_(1 : st, st + 1 : n)' ...
                            * B(1 : st, : )));
                        
                        %implement operatio A_'*B'
                    case 'T'
                        if(st > size(B, 2))
                            error('MESS:error_arguments', 'number of rows of A_ differs with cols of B');
                        end
                        C = (eqn.A_(1 : st, 1 : st)' + p * eqn.E_(1 : st, 1 : st)') * B( : , 1 : st)' ...
                            - eqn.A_(st + 1 : n, 1 : st)' ...
                            * (eqn.A_(st + 1 : n, st + 1 : n)' \ (eqn.A_(1 : st, st + 1 : n)' ...
                            * B( : , 1 : st)'));
                end
        end
        
end
end
