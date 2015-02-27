function [ mse ] = psnr( a,b )
%PSNR Peak-signal-to-noise ratio

if size(a,3) > 1 || size(b,3) > 1,
    warning('psnr only on grayscale images');
    mse = 0;
    return;
end

[m,n] = size(a);
if m ~= size(b,1) || n ~= size(b,2),
    error('images do not have same size');
end

ad = double(a);
bd = double(b);

err = ad-bd;
mse = sum(sum(err .* err ));
mse = mse / (m*n);

if mse > 0,
    mse = 10*log(255*255 / mse) / log(10);
else
    mse = 99;
end


end

