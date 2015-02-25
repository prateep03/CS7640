function [ fx,fy ] = gradient( F, options )
% GRADIENT - gradient operator, forward differences
%
% [fx,fy] = gradient(F, options)
% or
% f = gradient(F,options)

options.empty = 0;

fx = F([2:end end],:) - F;
fy = F(:,[2:end end]) - F;

if nargout == 1,
    fx = cat(3,fx,fy);
end

end

