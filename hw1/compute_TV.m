function TV = compute_TV(X, options)

options.empty = 0;

ndim = 2;
if size(X,1) == 1 || size(X,2) == 1,
    ndim = 1;
end
if size(X,1) > 1 && size(X,2) > 1 && size(X,3) > 1,
    ndim = 3;
end

if ndim == 1, % 1D
    TV = sum(abs(diff(X)));
    return;
end

tv_norm = readParam(options,'tv_norm','l2');

g = gradient(X,options);

switch tv_norm,
    case 'l1'
        TV = sum( abs(g),ndim+1 );
    case 'l2'
        TV = sqrt( sum(g.^2,ndim+1) );
    otherwise
        error('Not supported')
end
TV = sum(TV(:));