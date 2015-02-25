function [ Xpm, errImg, err ] = denoise_PM( X, options )
%DENOISE_PM Denoise image with Perona-Malik isotropic diffusion.
%   Detailed explanation goes here
% VARARGIN: Optional parameters:
%       sigma : standard deviation value
%       time : Time parameter 
%       nIter : Maximum number of iterations

% P. Perona, J. Malik:	Scale-Space and Edge Detection Using Anisotropic
% Diffusion, IEEE Trans. Pattern Anal. Mach. Intell., 
% 12(7), 1990, 629?639.

options.empty = 0;

ndim = getDim(X);

sigma = readParam(options,'sigma',20);
time = readParam(options,'time',3);
niter = readParam(options,'niter',1000);
prepr = readParam(options,'preSmooth',0);
lambda = readParam(options,'lambda',0.25);
gaussSigma = readParam(options,'gaussSigma',0.5);
disp = readParam(options,'display',1);

if or(lambda < 0, lambda > 2/8),
    warning('lambda too high or negative.');
end

% Adding 1-border to image to avoid boundary problems
X = [X(:,1), X, X(:,size(X,2))];
X = vertcat(X(1,:), X, X(size(X,1),:));

% If initial solution is provided, use it, otherwise use default image
Xpm = readParam(options,'initialImage',[]);
if isempty(Xpm),
   Xpm = X;
else
   Xpm = [Xpm(:,1), Xpm, Xpm(:,size(Xpm,2))];
   Xpm = vertcat(Xpm(1,:), Xpm, Xpm(size(Xpm,1),:));
end

% stencil initialization
stencilN = zeros(size(X,1),size(X,2));
stencilS = zeros(size(X,1),size(X,2));
stencilE = zeros(size(X,1),size(X,2));
stencilW = zeros(size(X,1),size(X,2));
stencilNorm = zeros(size(X,1),size(X,2));

% Residual
errImg = X; %zeros(size(X,1),size(X,2));
err = 1e7;
errarr = [1e7 1e7 1e7];

%% Pre-smooth image for ease of gradient computation 
if ~prepr,
    h = fspecial('gaussian', round(2*gaussSigma+1), gaussSigma);
end

iter = 0;
ndisp = max(round(niter/100),2);
while ((sum(errarr) > 0) && (iter < niter)),
    % smooth image
    if ~prepr,
        Xpm = imfilter(Xpm,  h, 'symmetric');
    end
    coeff = edgeStop1(Xpm, sigma,options);
    coeff = coeff * time;
    
    % Stencil computation
    stencilN(2:end-1,2:end-1) = (coeff(2:end-1, 2:end-1) + coeff(1:end-2, 2:end-1))/2;
    stencilS(2:end-1,2:end-1) = (coeff(2:end-1, 2:end-1) + coeff(3:end, 2:end-1))/2;
    stencilE(2:end-1,2:end-1) = (coeff(2:end-1, 2:end-1) + coeff(2:end-1, 3:end))/2;
    stencilW(2:end-1,2:end-1) = (coeff(2:end-1, 2:end-1) + coeff(2:end-1, 1:end-2))/2;
    
    stencilNorm = stencilN + stencilS + stencilW + stencilE + 1;
    
    % compute solution using Gauss Seidel
    Xpm(2:2:end-1,2:2:end-1) = (X(2:2:end-1, 2:2:end-1) ...
        + lambda*(stencilN(2:2:end-1,2:2:end-1) .* Xpm(1:2:end-2,2:2:end-1) ...
        +  stencilS(2:2:end-1,2:2:end-1) .* Xpm(3:2:end,2:2:end-1) ...
        +  stencilE(2:2:end-1,2:2:end-1) .* Xpm(2:2:end-1,3:2:end) ...
        +  stencilW(2:2:end-1,2:2:end-1) .* Xpm(2:2:end-1,1:2:end-2))) ...
        ./ stencilNorm(2:2:end-1,2:2:end-1);
    
    Xpm(3:2:end, 3:2:end) = (X(3:2:end, 3:2:end) ...
        + lambda*(stencilN(3:2:end, 3:2:end) .* Xpm(2:2:end-1, 3:2:end) ...
        + stencilS(3:2:end, 3:2:end) .* Xpm(4:2:end, 3:2:end) ...
        + stencilE(3:2:end, 3:2:end) .* Xpm(3:2:end, 4:2:end) ...
        + stencilW(3:2:end, 3:2:end) .* Xpm(3:2:end, 2:2:end-1))) ...
        ./ stencilNorm(3:2:end, 3:2:end);

    Xpm(2:2:end-1, 3:2:end) = (X(2:2:end-1, 3:2:end) ...
        + lambda*(stencilN(2:2:end-1, 3:2:end) .* Xpm(1:2:end-2, 3:2:end) ...
        + stencilS(2:2:end-1, 3:2:end) .* Xpm(3:2:end, 3:2:end) ...
        + stencilE(2:2:end-1, 3:2:end) .* Xpm(2:2:end-1, 4:2:end) ...
        + stencilW(2:2:end-1, 3:2:end) .* Xpm(2:2:end-1, 2:2:end-1))) ...
        ./ stencilNorm(2:2:end-1, 3:2:end);

    Xpm(3:2:end, 2:2:end-1) = (X(3:2:end, 2:2:end-1) ...
        + lambda*(stencilN(3:2:end, 2:2:end-1) .* Xpm(2:2:end-1, 2:2:end-1) ...
        + stencilS(3:2:end, 2:2:end-1) .* Xpm(4:2:end, 2:2:end-1) ...
        + stencilE(3:2:end, 2:2:end-1) .* Xpm(3:2:end, 3:2:end) ...
        + stencilW(3:2:end, 2:2:end-1) .* Xpm(3:2:end, 1:2:end-2))) ...
        ./ stencilNorm(3:2:end, 2:2:end-1);
       
    % compute residual 
    errImg(2:end-1, 2:end-1) = ...
        -( stencilN(2:end-1,2:end-1) .* Xpm(1:end-2,2:end-1) ...
        +  stencilS(2:end-1,2:end-1) .* Xpm(3:end, 2:end-1) ...
        +  stencilE(2:end-1,2:end-1) .* Xpm(2:end-1, 3:end) ...
        +  stencilW(2:end-1,2:end-1) .* Xpm(2:end-1, 1:end-2)) ...
        +  stencilNorm(2:end-1,2:end-1) .* Xpm(2:end-1,2:end-1) - X(2:end-1,2:end-1);
    
    errnew = sum( sum(real(errImg).^2) );
    errdiff = errnew - err;
    err = errnew;
    errarr = [errdiff errarr(1, 1:(size(errarr,2)-1))];
    
    % Duplicate edges 
    Xpm = [Xpm(:,2), Xpm(:,2:end-1), Xpm(:,end-1)];
    Xpm = vertcat(Xpm(2,:), Xpm(2:end-1,:), Xpm(end-1,:));
    
    if mod(iter,ndisp) == 0 && ~isempty(disp),
         if ndim == 1,
             plot(Xpm); axis tight;
         else
             imageplot(Xpm);
         end
         drawnow;
     end
    
    iter = iter + 1;
    fprintf('Iteration : %d, residual : %.10f\n',iter,errdiff);
end

% Remove borders
Xpm = Xpm(2:end-1,2:end-1);
if nargout > 1,
    errImg = errImg(2:(size(Xpm,1)-1), 2:(size(Xpm,2)-1));
end


function d = getDim(x)
% getDim(x) : get dimensions of x(1D / 2D)
if isempty(x),
    d = 0;
    return;
end
d = ndims(x);
if d==2 && (size(x,1) == 1 || size(x,2) == 1)
    d = 1;
end

function f = edgeStop(im, sigma,options)
% PERONAEDGESTOP Perona Edge-Stopping function on matrix
% img: 2D image matrix

[fx,fy] = gradient(im,options);
f = sqrt(fx.*fx + fy.*fy);
f = exp(- f .* f / (2 * sigma * sigma));

function f = edgeStop1(im, sigma, options)
% img: 2D image matrix

[fx,fy] = gradient(im,options);
f = sqrt(fx.*fx + fy.*fy);
f = 1 ./ (1 + (f ./ sigma).^2);