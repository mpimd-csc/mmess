function [ y, perm ] = mess_make_proper( x )
%MESS_MAKE_PROPER ensures the input vector to be a proper set of shifts
%
%  The function checks if the vector contains set of shifts closed under
%  complex conjugation and all complex pairs are subsequent entries.
if isempty(x), y=x; return; end

[m,n] = size(x);

if m < n, x = x'; end

idro = find(imag(x)==0.0);
idco = find(imag(x));

xr = x(idro);
xc = x(idco);

k = length(xc);
if rem(k,2)~=0, error('odd number of complex shifts detected'); end

% Sort real shifts
[xr,idr] = sort(xr);

% Sort complex shifts w.r.t. real part
[~,idcr] = sort(real(xc));
xc = xc(idcr); 
% Complex conjugated pairs ensured
xc(2:2:k,:) = conj(xc(1:2:k,:));    

% % Sort complex shifts w.r.t. complex part
% [~, idcc]=sort(imag(xc),'descend');
% xc = xc(idcc); 

% Sort complex shifts w.r.t. real part
[~, idcr2]=sort(real(xc));
xc = xc(idcr2);

y = [xr; xc];
perm = [idro(idr);idco(idcr(idcr2))];
end

