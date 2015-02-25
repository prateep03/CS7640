M0 = load_image('img1_orig',256);
if(size(M0,3) > 1),
    tmp = zeros(size(M0,1),size(M0,2));
    tmp = (M0(:,:,1) + M0(:,:,2) + M0(:,:,3) ) / 3;
    M0 = tmp;
end
sigma = .1;
n = 256;
M0 = rescale(crop(M0,n));
M = load_image('img1_sigma_25',256); 
% M = M0 + randn(size(M0))*sigma;
if(size(M,3) > 1),
    tmp = zeros(size(M0,1),size(M0,2));
    tmp = (M(:,:,1) + M(:,:,2) + M(:,:,3) ) / 3;
    M = tmp;
end
M = rescale(crop(M,n));
options.verb = 0;
options.display = 0;
options.niter = 50;    % number of iterations
options.niter_inner = 100;
options.lambda = .1; % initial regularization

[Mtv,err,tv] = chambolle(M,options);
imageplot({Mtv M0 M},{'denoised','original','noisy'});

figure; clf;
lagr = .5*err.^2 + options.lambda*tv;
h = plot(lagr); axis tight;
set(h, 'LineWidth', 2);
set(gca, 'FontSize', 20);
title('Evolution of |y-Ux|^2+\lambda |x|_1');
axis([1 options.niter min(lagr) min(lagr)*1.3]);