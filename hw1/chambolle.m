function [ Xtv, err, tv ] = chambolle( x, options )
% Chambolle denoising
%       A. Chambolle: An Algorithm for Total Variation Minimization and Applications.
%       Journal of Mathematical Imaging and Vision 20 (1-2): 89-97, January - March, 2004

options.empty = 0;

ndim = getDim(x);

%% Get all parameters
la = readParam(options, 'lambda', []);
if isempty(la),
    la = .1;
end

verbose = readParam(options,'verbose',1);
niter = readParam(options,'niter',3000); 
niter_inner = readParam(options,'niter_inner',30);
tol = readParam(options,'tol',1e-5);
disp = readParam(options,'display',1);

% "In practice, it appears that the optimal constant for the stability
% and convergence of the algorithm is not 1/8 but 1/4".
tau = readParam(options,'tau',1/4);
if tau > 1/4,
    warning('tau is too large for convergence');
    tau = 1/4;
end

x1 = x; % initialize
% initialize p
if ndim == 1, % 1D
    p = zeros(size(x));
else % 2D
    p = zeros([size(x) 2]);
end
p = readParam(options,'p',p);
TV = readParam(options,'TV', []);

if isempty(TV),
    niter_inner = 1;
end

if disp,
    clf;
end

err = []; tv = [];
ndisp = max(round(niter/100),2);
for i=1:niter
     fprintf('Iteration : %d\n', i);    
     for ii=1:niter_inner
         % chambolle step
         gp = gradient( divergence(p,options) - x/la, options );
         if ndim == 1,
             d = abs(g);
         else
             d = sqrt(sum(gp.^2,3));
             d = repmat(d, [1 1 2]);
         end
         pnew = (p + tau*gp) ./ (1 + tau*d);
         x2 = x - la*divergence(pnew,options);
         rel_err = norm(x2 - x1,'fro') / sqrt(numel(x1));
         if rel_err < tol, % converged
             x1 = x2;
             break;
         end
         x1 = x2; % update solution
         p = pnew;
     end
     
     if nargout >= 2,
         if isempty(TV),
             err(i) = norm(x - x1, 'fro');  
         end
     end
     if nargout >= 3,
         tv(i) = compute_TV(x1, options);
     end
          
     if mod(i,ndisp) == 0 && ~isempty(disp),
         if ndim == 1,
             plot(x1); axis tight;
         else
             imageplot(x1);
         end
         drawnow;
     end
end
Xtv = x1;

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