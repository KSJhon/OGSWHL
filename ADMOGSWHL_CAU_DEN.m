function out = ADMOGSWHL_CAU_DEN( Bn, r, params)
%beta1, beta2, maxItr, K, p, sigma_w
[m, n] = size(Bn);
resfg  = 1.e-5;
K1      = params.K;
K2      = params.K;
kk      = 0;
NIt     = 20;
kkmax= 5;
tempfun = 0;
maxiter = params.MaxIter;
lambda = params.lambda;
q = params.q;
q2 = q * 2;
beta1 = params.beta(1);beta2 = params.beta(2);

f  = min(max(Bn, 0), 255); 
u1 = zeros(m, n);
w  = u1;
u2 = u1;
u3 = u1;

Demmo = beta1 * abs(psf2otf([1, -1], [m, n])).^2 + beta1 * abs(psf2otf([1; -1], [m, n])).^2 + beta2;  
otf1  = psf2otf([1, -1], [m, n]);
otf2  = psf2otf([1; -1], [m, n]);

% finite diff
[D, ~] = defDDt;
[DuX, DuY] = D(f);
h = fspecial('gaussian', [K1 K2], params.sigma_w);
s_w = h * K1 * K2;
s_w2 = s_w.^2;

gs1 = sqrt(conv2(abs(DuX).^q2, s_w2, 'same'));
gs2 = sqrt(conv2(abs(DuY).^q2, s_w2, 'same'));

[g] = fval(gs1, gs2, f, Bn, lambda, r);

out.g = g;
inv_beta1 = 1/beta1;
%% Main loop
for ii = 1:maxiter
    
    % ==================
    %     v-subprolem
    % ==================
    d1 = DuX + u1 * inv_beta1;
    d2 = DuY + u2 * inv_beta1;
    
    v1 = gstvdm_weighted(d1, K1, K2, inv_beta1, NIt, q, s_w2);
    v2 = gstvdm_weighted(d2, K1, K2, inv_beta1, NIt, q, s_w2);

    % ==================
    %     w-subprolem
    %     Newton method
    % ==================
    temp = r^2 + (w - Bn).^2;
    fd   = lambda * ((w - Bn)./temp) + beta2 * (w - f) - u3;
    while kk < kkmax && (norm(fd, 'fro') > 1e-8)
        fdd = lambda * ((r^2 - (w - Bn).^2)./(temp.^2)) + beta2;
        w   = w - fd./fdd;
        kk  = kk + 1;
        
        temp = r^2 + (w - Bn).^2;
        fd   = lambda*((w - Bn)./temp) + beta2 * (w - f) - u3;
    end
    kk = 0;
    w(w < 0) = 0;
    
    % ==================
    %     f-subprolem
    % ==================
    nomin1 = conj(otf1).*fft2(v1) + conj(otf2).*fft2(v2);
    nomin2 = conj(otf1).*fft2(u1) + conj(otf2).*fft2(u2);
    
    FW = beta1 * nomin1 + beta2 * fft2(w) - nomin2 - fft2(u3);
    
    f = real(ifft2(FW./Demmo));
    
    % finite diff.
    [DuX, DuY] = D(f);
    
    gs1 = sqrt(conv2(abs(DuX).^q2, s_w2, 'same'));
    gs2 = sqrt(conv2(abs(DuY).^q2, s_w2, 'same'));
    [g] = fval(gs1, gs2, f, Bn, lambda, r);
    
    out.g = [out.g; g];
    % ====================
    % Check stopping rule
    % ====================

    if abs(out.g(ii) - tempfun)/abs(tempfun) < resfg
        out.sol = f;
        [DuX,DuY] = D(f);
        
        gs1 = sqrt(conv2(abs(DuX).^q2, s_w2, 'same'));
        gs2 = sqrt(conv2(abs(DuY).^q2, s_w2, 'same'));
        [g] = fval(gs1, gs2, f, Bn, lambda, r);
        out.g = [out.g; g];
        return
    end
    
    tempfun = out.g(ii);
    % ==================
    %    Update Lam
    % ==================
    u1 = u1 + beta1 * (DuX - v1);
    u2 = u2 + beta1 * (DuY - v2);
    u3 = u3 + beta2 * (f - w);
end

out.sol = f;
if ii == maxiter
    out.exit = 'Maximum iteration reached!';
end

%% ------------------SUBFUNCTION-----------------------------
function [g] = fval(gs1, gs2, f, Bn, lambda, r)

reg = sum(sum(gs1 + gs2));

fid = (lambda/2) * sum(sum(log(r^2 + (f - Bn).^2)));
g = reg +  fid;

function [D, Dt] = defDDt

D = @(U) ForwardD(U);
Dt = @(X, Y) Dive(X, Y);

function [Dux, Duy] = ForwardD(U)

Dux = [diff(U, 1, 2), U(:, 1) - U(:, end)];
Duy = [diff(U, 1, 1); U(1, :) - U(end, :)];

function DtXY = Dive(X, Y)

DtXY = [X(:, end) - X(:, 1), -diff(X, 1, 2)];
DtXY = DtXY + [Y(end, :) - Y(1, :); -diff(Y, 1, 1)];
