function x  = gstvdm_weighted(y, K1, K2, lam, Nit, p, s_w2)

x = y;                                              % Initialization

p2 = p * 2;
if K1 ~= 1 || K2 ~= 1
    for k = 1:Nit
        r = sqrt(conv2(abs(x).^p2, s_w2, 'same')) + eps; % zero outside the bounds of  x
        v = conv2(1./r, s_w2, 'same');
        x = y./(1 + p * lam * v.*abs(x).^(p2 - 2));
    end
else
    for k = 1:Nit
        r = sqrt(abs(x).^p2); 
        v = 1./r;  
        x = y./(1 + lam * v);
     end
end

