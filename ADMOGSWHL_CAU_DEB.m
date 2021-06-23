function out = ADMOGSWHL_CAU_DEB(H, Bn, r, params)
%.lambda, mu, delta1, delta3, delta2, K, p, sigma_h, maxItr
[m, n] = size(Bn);
I0 = medfilt2(min(max(Bn, 0), 1));

resfg = 1.e-5;

% initialization
maxiter = params.MaxIter;
K1 = params.K;
K2 = params.K;
lambda = params.lambda;
mu = params.mu;

NIt = 20;
kkmax = 5;
kk = 0;
tempfun = 0;

f = Bn;
w1 = zeros(m, n);
y = w1;
w2 = w1;
w4 = w1;
w3 = w1;

otf = psf2otf(H, [m, n]);
otf1 = psf2otf([1, -1], [m, n]);
otf2 = psf2otf([1; -1], [m, n]);
delta1 = params.beta(1);delta2 = params.beta(3);delta3 = params.beta(2);
Demmo = delta1 * abs(otf1).^2 + delta1 * abs(otf2).^2 + delta3 * abs(otf).^2 + delta2;


% finite diff
[D, ~] = defDDt;
[DuX, DuY] = D(f);

q = params.q;
q2 = q * 2;
h = fspecial('gaussian', [K1 K2], params.sigma_w);
s_h = h * K1 * K2;
s_h2 = s_h.^2;

gs1 = sqrt(conv2(abs(DuX).^q2, s_h2, 'same'));
gs2 = sqrt(conv2(abs(DuY).^q2, s_h2, 'same'));

KF = imfilter(f, H, 'circular', 'conv');

[g] = fval(gs1, gs2, KF, Bn, lambda, mu, I0, r);

out.g = g;
%% Main loop
for ii = 1:maxiter
    % ==================
    %     x1, x2-subprolem
    % ==================
    d1 = DuX + w1 / delta1;
    d2 = DuY + w2 / delta1;
    
    v1 = gstvdm_weighted(d1, K1, K2, 1 / delta1, NIt, q, s_h2);
    v2 = gstvdm_weighted(d2, K1, K2, 1 / delta1, NIt, q, s_h2);
 
     % ==================
    %     y-subprolem
    %     Newton method
    % ==================
    temp = r^2 + (y - Bn).^2;
    fd = lambda * ((y - Bn)./temp + mu.*(y - I0)) + delta3 * (y - KF) - w4;
    while kk < kkmax && (norm(fd, 'fro')>1e-8)
        
        fdd = lambda * ((r^2 - (y - Bn).^2)./(temp.^2) + mu) + delta3;
        y = y - fd./fdd;
        kk = kk + 1;
        
        temp = r^2 + (y - Bn).^2;
        fd = lambda * ((y - Bn)./temp + mu.*(y - I0)) + delta3 * (y - KF) - w4;
    end
    kk = 0;
    y(y < 0) = 0;
    %===================
    %    x3-subproblem
    %===================
    
    x3 = min(1, max(f + w3/delta2, 0));
    
    % ==================
    %     f-subprolem
    % ==================
   
    nomin1 = conj(otf1).*( delta1 * fft2(v1) - fft2(w1)) + conj(otf2).*( delta1 * fft2(v2) - fft2(w2));
    nomin2 = conj(otf).*( delta3  * fft2(y) - fft2(w4));
    
    FW = nomin1 + nomin2 + delta2 * fft2(x3) - fft2(w3);
    
    f = real(ifft2(FW./Demmo));
    
    [DuX, DuY] = D(f);
    
    gs1 = sqrt(conv2(abs(DuX).^q2, s_h2, 'same'));
    gs2 = sqrt(conv2(abs(DuY).^q2, s_h2, 'same'));
    KF = imfilter(f, H, 'circular', 'conv');
    
    [g] = fval(gs1, gs2, KF, Bn, lambda, mu, I0, r);
    out.g = [out.g; g];

    if abs(out.g(ii) - tempfun) / abs(tempfun) < resfg
        out.sol = f;
        out.itr = ii;
        [DuX, DuY] = D(f);
        
        gs1 = sqrt(conv2(abs(DuX).^q2, s_h2, 'same'));
        gs2 = sqrt(conv2(abs(DuY).^q2, s_h2, 'same'));
        KF = imfilter(f, H, 'circular', 'conv');
        
        [g] = fval(gs1, gs2, KF, Bn, lambda, mu, I0, r);
        out.g = [out.g; g];
       
        return
    end
    
    tempfun = out.g(ii);
    % ==================
    %    Update Lam
    % ==================
    w1 = w1 + delta1 * (DuX - v1);
    w2 = w2 + delta1 * (DuY - v2);
    w3 = w3 + delta2 * (f - x3);
    w4 = w4 + delta3 * (KF - y);
end

out.sol = f;
out.itr = ii;
out.exit = 'Exist Normally';
if ii == maxiter
    out.exit = 'Maximum iteration reached!';
end


%% ------------------SUBFUNCTION-----------------------------
function [g] = fval(gs1, gs2, KF, Bn, lambda, mu, I0, r)

tv =  sum(sum(gs1 +gs2));

fid =  (lambda / 2) * sum(sum(log(r^2 + (KF - Bn).^2) + mu * (KF - I0).^2));
g = tv +  fid;

function [D, Dt] = defDDt

D = @(U) ForwardD(U);
Dt = @(X, Y) Dive(X, Y);

function [Dux, Duy] = ForwardD(U)

Dux = [diff(U, 1, 2), U(:, 1) - U(:, end)];
Duy = [diff(U, 1, 1); U(1, :) - U(end, :)];

function DtXY = Dive(X, Y)

DtXY = [X(:, end) - X(:, 1), -diff(X, 1, 2)];
DtXY = DtXY + [Y(end, :) - Y(1, :); -diff(Y, 1, 1)];

function x = randl(k, sigma)
for i = 1:k
    for j = i:k
        x(i, j) = exp(-(abs(i - ceil(k / 2)) + abs(j - ceil(k/2))) / sigma);
    end
end
for i = 1:k
    for j = 1:i - 1
        x(i, j) = x(j, i);
    end
end
x = x./sum(x(:));