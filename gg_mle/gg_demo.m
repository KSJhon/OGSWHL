% GGDEMO 
% Demonstrations of the main functions in the generalized Gaussian density package 
% 
 figure
 
disp('Model parameter:') 
mu = 0       % mean 
alpha = 4     % scale 
beta = .8    % shape 
 
% Generate 10^4 random samples from the generalize Gaussian density 
r = ggrnd(mu, alpha, beta, 1, 8*10^4); 
 
disp('Moment matching estimate:'); 
[mu1, alpha1, beta1] = ggmme(r) ;
 
disp('Maximum likelihood estimate:'); 
[mu2, alpha2, beta2] = ggmle(r) ;
         disp(sprintf('[mu2=%.3f,alpha2=%.3f,beta2=%.3f]',mu2, alpha2, beta2));

[mu22, alpha22, beta22] = ggmle_modified(r) ;
        disp(sprintf('[mu2=%.3f,alpha2=%.3f,beta2=%.3f]',mu22, alpha22, beta22));

% Compare the estimated PDF's with the histogram 
[N, X] = hist(r, 31); 
 
clf; 
bar(X, N ./ (X(2) - X(1)) / sum(N)); 
hold; 
 
plot(X, ggpdf(X, mu1, alpha1, beta1), 'b'); 
plot(X, ggpdf(X, mu2, alpha2, beta2), 'r'); 
plot(X, ggpdf(X, mu22, alpha22, beta22), 'g'); 