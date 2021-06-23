% Weighted Hyper-Laplacian Prior with Overlapping Group Sparsity for Image Restoration under Cauchy Noise
% written by K.S. Jon, 20200426
% Matlab Version 9.1.0.441655 (R2016b)
%% Run OGS-WeightedHL deblurring
clear variables
clc
close all

fileList = {'baby.tif', 'boat.pgm', 'lena.pgm', 'parrot.png', 'camman.png', 'WStatn.tif', 'jelfih.png', 'wo_dhr.tif'};

% some constants
r  = 0.02; % noise level
params.q = 0.8; % exponent for hyper-Laplacian
params.sigma_w = 4; % weighting parameter
params.K  = 5; % group size for OGS
params.MaxIter = 500; % maimum iteration number 
params.beta = [420 280 60]; % penalty parameter for ADMM

if 0 	% Gaussian deblurring simulation
    H = fspecial('gaussian', [7 7], 3);
    mu_list = [19    24    24    24    24    26    28    14];
    
    lamda_list = [37    56    38    61    63    52    33    23];
else	% motion deblurring simulation
    H = fspecial('motion', 8, 30);
    mu_list = [8    22     8     8     7    20    20    16];
    lamda_list = [32    41    33    40    46    35    25    18];
end
for j =1:size(fileList, 2)
    img_file = fileList{j};
    I = double(imread(img_file))./255;
   
    % first, mimic blurring effect
    Bn = imfilter(I, H, 'circular', 'conv');
    % second, add Cauchy noise
    rng('default')
    v1 = randn(size(I));
    v2 = randn(size(I));
    Bn = Bn + r * v1./v2;
    Bn = min(max(Bn, -10), 20); % preprocessing
    % calculate quality metric for blurry image
    psnr_blurry = psnr(I, min(max(Bn, 0), 1), 1);
    ssim_blurry = ssim(min(max(Bn, 0), 1) * 255, I * 255);
    
    params.mu = mu_list(j);
    params.lambda = lamda_list(j);
    tic;
    [out] = ADMOGSWHL_CAU_DEB(H, Bn, sqrt(r), params);
    time = toc;
    x  = min(max(out.sol, 0), 1);
    % calculate quality metric for recovered image
    psnr_recon  = psnr(x, I, 1);
    ssim_recon = ssim(x * 255, I * 255);
    
    display(sprintf('%s:\tpsnr0=%1.2f\tssim0=%1.3f\tpsnr1=%1.2f\tssim1=%1.3f\ttime=%.2f\n', ...
        img_file, psnr_blurry, ssim_blurry, psnr_recon,  ssim_recon, time));
	display(sprintf('=================================='));
    
    figure;
    subplot(1, 2, 1); imshow(min( max( Bn, 0 ), 1), []); title({'Noisy image';['PSNR = ', num2str(psnr_blurry, '%3f')];['SSIM = ', num2str(ssim_blurry, '%2f')]});
    subplot(1, 2, 2); imshow(x, []); title({'Recovered image';['PSNR = ', num2str(psnr_recon, '%3f')];['SSIM = ', num2str(ssim_recon, '%2f')]});
end