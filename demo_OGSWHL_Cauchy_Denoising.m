% Weighted Hyper-Laplacian Prior with Overlapping Group Sparsity for Image Restoration under Cauchy Noise
% written by K.S. Jon, 20200426
% Matlab Version 9.1.0.441655 (R2016b)
%% Run weighted OGSHL for Cauchy noise removal
clear variables;
close all;
clc;
gamma = 5; 	% noise level
% some constants
params.K  = 5;		% group size for OGS
params.q = .8;		% exponent for hyper-Laplacian
params.beta = [0.3 22]; % penalty parameter for ADMM
params.sigma_w = 4;	% weighting parameter
params.MaxIter = 500;	% maimum iteration number 

fileList = {'baby.tif', 'boat.pgm', 'lena.pgm', 'parrot.png', 'cama(256).png', 'WStatn.tif', 'jelfih.png', 'wo_dhr.tif'};
if gamma == 5
    lambdaList = [84    95    90    94    97    90    83    80];
else %gamma == 10
    lambdaList = [134   146   137   140   145   139   138   131];
end
display(sprintf('OGSWHL_Denoise_%d', gamma));

for j =1:size(fileList,2)
    img_file = fileList{j};
   
    I = double(imread(img_file));
    %% add Cauchy noise to clean image

    randn('state', -34);
    v1 = randn(size(I));
    randn('state', 94);
    v2 = randn(size(I));
    Bn = I + gamma * v1./v2;

    psnr_noisy = psnr(min(max(Bn, 0), 255), I, 255);
    ssim_noisy = ssim(min(max(Bn, 0), 255), I);

    params.lambda = lambdaList(j);
    tic;
	[out] = ADMOGSWHL_CAU_DEN( Bn, gamma, params);
    time = toc;

	x = min(max(out.sol,0), 255);
	psnr_recon = psnr(x, I, 255);
	[ssim_recon ssim_map]= ssim(x, I);
	
	display(sprintf('%s:\tpsnr0=%1.2f\tssim0=%1.3f\tpsnr1=%1.2f\tssim1=%1.3f\ttime=%.2f\n', ...
        img_file, psnr_noisy, ssim_noisy, psnr_recon,  ssim_recon, time));
	display(sprintf('=================================='))
	
	figure;
    subplot(1, 2, 1); imshow(min(max(Bn, 0), 255), []); title({'Noisy image';['PSNR = ', num2str(psnr_noisy, '%3f')];['SSIM = ', num2str(ssim_noisy, '%2f')]});
    subplot(1, 2, 2); imshow(x, []); title({'Recovered image';['PSNR = ', num2str(psnr_recon, '%3f')];['SSIM = ', num2str(ssim_recon, '%2f')]});
end
