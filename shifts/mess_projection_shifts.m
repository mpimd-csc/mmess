function p = mess_projection_shifts(eqn, opts, oper, V, W, p_old)
%%function p = mess_projection_shifts(eqn, opts, oper, V, W, p_old)
%
% Internal helper function for usfs and mess_get_projection
% shifts. Computes new shifts by implicitly or explicitly
% projecting the E and A matrices to the span of V. Note that the
% width of V must be a multiple of that of W, V is the newest part
% of the ADI solution factor Z and the old shift
% vector p_old passed in must have this multiple as its length.
%
% Whether or not the projection is computed implicitly from the
% contents of V or by an explicit projection, is determined via
% opts.shifts.implicitVtAV.
%

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%


%% Check data
if not(isfield(opts,'shifts')) || not(isstruct(opts.shifts))
    warning('MESS:control_data',['shift parameter control structure missing.', ...
        'Switching to default num_desired = 25.']);
    opts.shifts.num_desired = 25;
else
    if not(isfield(opts.shifts,'num_desired')) || ...
       not(isnumeric(opts.shifts.num_desired))
        warning('MESS:control_data',...
            ['Missing or Corrupted opts.shifts.num_desired field.', ...
            'Switching to default: 25']);
        opts.shifts.num_desired = 25;
    end
    if not(isfield(opts.shifts,'implicitVtAV'))|| isempty(opts.shifts.implicitVtAV)
        opts.shifts.implicitVtAV = true;
    end
end

if not(isfield(eqn, 'haveE')), eqn.haveE = 0; end
[result, eqn, opts, oper] = oper.init(eqn, opts, oper, 'A', 'E');
if not(result)
    error('MESS:control_data', 'system data is not completely defined or corrupted');
end

L = length(p_old);
nV = size(V, 2);
nW = size(W, 2);
if L > 0 && any(p_old)
    if nV / nW ~= L
        error('MESS:control_data', 'V and W have inconsistent no. of columns');
    end
end

%% Initialize data
if L > 0 && any(p_old)
    T = zeros(L, L);
    K = zeros(1, L);
    D = [];
    Ir = eye(nW);
    iC = find(imag(p_old));
    iCh = iC(1 : 2 : end);
    iR = find(not(imag(p_old)));
    isubdiag = [iR; iCh];
    h = 1;
end

%% Process previous shifts
if L > 0 && any(p_old) && opts.shifts.implicitVtAV
    while h <= L
        is = isubdiag(isubdiag < h);
        K(1, h) = 1;
        if isreal(p_old(h)) % real shift
            T(h, h) = p_old(h);
            if not(isempty(is))
                T(h, is) = 2 * p_old(h) * ones(1, length(is));
            end
            D = blkdiag(D, sqrt(-2 * p_old(h)));
            h = h + 1;
        else % complex conjugated pair of shifts
            rpc=real(p_old(h));
            ipc=imag(p_old(h));
            beta=rpc / ipc;
            T(h : h + 1, h : h + 1) = [3 * rpc, -ipc;
                                       ipc * (1 + 4 * beta^2), -rpc];
            if not(isempty(is))
                T(h : h+  1, is)=[4 * rpc;
                                  4 * rpc * beta] * ones(1, length(is));
            end
            D = blkdiag(D, 2 * sqrt(-rpc) * [1,0; beta, sqrt(1 + beta^2)]);
            h = h + 2;
        end
    end
    S = kron(D \ (T * D), Ir);
    K = kron(K * D, Ir);
else  % explicit AV (unless already computed in mess_para)
    S = 0;
    K = 1;
    if any(p_old)
        W = oper.mul_A(eqn, opts, eqn.type, V, 'N');
        if isfield(eqn, 'haveUV') && eqn.haveUV
            if eqn.type == 'T'
                W = W + eqn.V * (eqn.U' * V);
            else
                W = W + eqn.U * (eqn.V' * V);
            end
       end
    end
end

%% Compute projection matrices
[~, s, v] = svd(V' * V);
s = diag(s);
r = sum(s > eps * s(1) * nV);
st = v( : , 1 : r) * diag(1 ./ s(1 : r).^.5);
U = V * st;


%% Project V and compute Ritz values
if eqn.haveE
    E_V = oper.mul_E(eqn, opts, eqn.type, V, 'N');
    G = U' * E_V;
    H = U' * W * K * st + G * (S * st);
    G = G * st;
    p = eig(H, G);
else
    H = U' * (W * K) * st + U'*( V *( S * st ));
    p = eig(H);
end

%% Postprocess new shifts

% remove infinite values
p = p(isfinite(p));
% remove zeros
p = p(abs(p) > eps);
% make all shifts stable
p(real(p) > 0) = -p(real(p) > 0);
if not(isempty(p))
    % remove small imaginary pertubations
    small_imag = find(abs(imag(p)) ./ abs(p) < 1e-12);
    p(small_imag) = real(p(small_imag));
    % sort (s.t. compl. pairs are together)
    p=sort(p);
    if length(p) > opts.shifts.num_desired
        p = mess_mnmx(p, opts.shifts.num_desired);
    end
end
