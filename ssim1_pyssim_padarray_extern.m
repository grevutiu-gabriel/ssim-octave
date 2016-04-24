function [mssim, ssim_map] = ssim1_pyssim_padarray_extern(img1_extern, img2, img1_squared_extern, mu11_extern, mu1_sq_extern, mu1111_extern)

% ========================================================================
% SSIM Index with automatic downsampling, Version 1.0
% Copyright(c) 2009 Zhou Wang
% All Rights Reserved.
%
% ----------------------------------------------------------------------
% Permission to use, copy, or modify this software and its documentation
% for educational and research purposes only and without fee is hereby
% granted, provided that this copyright notice and the original authors'
% names appear on all copies and supporting documentation. This program
% shall not be used, rewritten, or adapted as the basis of a commercial
% software or hardware product without first obtaining permission of the
% authors. The authors make no representations about the suitability of
% this software for any purpose. It is provided "as is" without express
% or implied warranty.
%----------------------------------------------------------------------
%
% This is an implementation of the algorithm for calculating the
% Structural SIMilarity (SSIM) index between two images
%
% Please refer to the following paper and the website with suggested usage
%
% Z. Wang, A. C. Bovik, H. R. Sheikh, and E. P. Simoncelli, "Image
% quality assessment: From error visibility to structural similarity,"
% IEEE Transactios on Image Processing, vol. 13, no. 4, pp. 600-612,
% Apr. 2004.
%
% http://www.ece.uwaterloo.ca/~z70wang/research/ssim/
%
% Note: This program is different from ssim_index.m, where no automatic
% downsampling is performed. (downsampling was done in the above paper
% and was described as suggested usage in the above website.)
%
% Kindly report any suggestions or corrections to zhouwang@ieee.org
%
%----------------------------------------------------------------------
%
%Input : (1) img1: the first image being compared
%        (2) img2: the second image being compared
%        (3) K: constants in the SSIM index formula (see the above
%            reference). defualt value: K = [0.01 0.03]
%        (4) window: local window for statistics (see the above
%            reference). default widnow is Gaussian given by
%            window = fspecial('gaussian', 11, 1.5);
%        (5) L: dynamic range of the images. default: L = 255
%
%Output: (1) mssim: the mean SSIM index value between 2 images.
%            If one of the images being compared is regarded as 
%            perfect quality, then mssim can be considered as the
%            quality measure of the other image.
%            If img1 = img2, then mssim = 1.
%        (2) ssim_map: the SSIM index map of the test image. The map
%            has a smaller size than the input images. The actual size
%            depends on the window size and the downsampling factor.
%
%Basic Usage:
%   Given 2 test images img1 and img2, whose dynamic range is 0-255
%
%   [mssim, ssim_map] = ssim(img1, img2);
%
%Advanced Usage:
%   User defined parameters. For example
%
%   K = [0.05 0.05];
%   window = ones(8);
%   L = 100;
%   [mssim, ssim_map] = ssim(img1, img2, K, window, L);
%
%Visualize the results:
%
%   mssim                        %Gives the mssim value
%   imshow(max(0, ssim_map).^4)  %Shows the SSIM index map
%========================================================================
format long
output_precision(7)

if (columns(size(img1_extern))!=2 || columns(size(img2))!=3)
   mssim = 123;
   ssim_map = 123;
   return;
end


%   window = fspecial('gaussian', 11, 1.5);	%
   gaussian_kernel_width=11;
   gaussian_kernel_sigma=1.5;
   window= rand(1,gaussian_kernel_width);
   norm_mu = floor(gaussian_kernel_width / 2);
   for i = 1:gaussian_kernel_width
   window(i)= ((1 / (sqrt(2 * pi) * (gaussian_kernel_sigma))) * exp(-((((i-1) - norm_mu) ** 2)) / (2 * (gaussian_kernel_sigma ** 2))));
   endfor
   K(1) = 0.01;					% default settings
   K(2) = 0.03;					%
   L = 255;
   %window                                     %

%img1 = rgb2gray(double(img1));
#img2 = rgb2gray(double(img2));
img1=img1_extern;
img2=floor(114/1000*double(img2(:,:,3))+587/1000*double(img2(:,:,2))+299/1000*double(img2(:,:,1))+0.00000001);

img1_squared = img1_squared_extern;
img2_squared = img2 .* img2;


C1 = (K(1)*L)^2;
C2 = (K(2)*L)^2;
% window = window/sum(sum(window));

mu11  = mu11_extern;
mu1_sq = mu1_sq_extern;

mu1111= mu1111_extern;
%mu1   = filter2(window, mu1, 'valid');
mu2   = filter2(window, padarray(img2', [0 5], "symmetric"))(:,6:1085)';
mu22  = filter2(window, padarray(mu2, [0 5], "symmetric"))(:,6:1925);
mu2_sq = mu22.*mu22;
mu222  = filter2(window, padarray(img2_squared', [0 5], "symmetric"))(:,6:1085)';
mu2222 = filter2(window, padarray(mu222, [0 5], "symmetric"))(:,6:1925);
%size(mu2)
%mu2   = filter2(window, mu2, 'valid');
sigma1_sq = mu1111;
sigma1_sq = sigma1_sq - mu1_sq;

sigma2_sq = mu2222;
sigma2_sq = sigma2_sq - mu2_sq;

img_mat_12 = img1 .* img2;
result12 = filter2(window, padarray(img_mat_12', [0 5], "symmetric"))(:,6:1085)';
result1212 = filter2(window, padarray(result12, [0 5], "symmetric"))(:,6:1925);
img_mat_sigma_12 = result1212;
img_mat_mu_12 = mu11.*mu22;
img_mat_sigma_12 = img_mat_sigma_12 - img_mat_mu_12;

%sigma1_sq = filter2(window, img1.*img1, 'valid') - mu1_sq;
%sigma2_sq = filter2(window, img2.*img2, 'valid') - mu2_sq;
%sigma12 = filter2(window, img1.*img2, 'valid') - mu1_mu2;

if (C1 > 0 && C2 > 0)
   ssim_map = ((2*img_mat_mu_12 + C1).*(2*img_mat_sigma_12 + C2))./((mu1_sq + mu2_sq + C1).*(sigma1_sq + sigma2_sq + C2));
%else
 %  numerator1 = 2*mu1_mu2 + C1;
 % numerator2 = 2*sigma12 + C2;
%	denominator1 = mu1_sq + mu2_sq + C1;
 %  denominator2 = sigma1_sq + sigma2_sq + C2;
  % ssim_map = ones(size(mu1));
   %index = (denominator1.*denominator2 > 0);
 %  ssim_map(index) = (numerator1(index).*numerator2(index))./(denominator1(index).*denominator2(index));
  % index = (denominator1 ~= 0) & (denominator2 == 0);
  % ssim_map(index) = numerator1(index)./denominator1(index);
end

mssim = mean2(ssim_map);

return
