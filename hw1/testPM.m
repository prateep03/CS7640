M = load_image('img1_sigma_12',256);
if(size(M,3) > 1),
    tmp = zeros(size(M,1),size(M,2));
    tmp = (M(:,:,1) + M(:,:,2) + M(:,:,3) ) / 3;
    M = tmp;
end

options.sigma = 20;
options.time = 3;
options.niter = 200;
options.preSmooth = 1;
options.display = 0;
options.lambda = 1.0;

[Mpm,errImg,err] = denoise_PM(M,options);
clf;
imageplot({Mpm M errImg},{'denoised','original','difference'});

tv = compute_TV(Mpm,options);
fprintf('Total variation = %.20e\n',tv);