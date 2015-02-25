function fd = divergence( Fx, Fy, options )

% backward difference

options.empty = 0;
if size(Fx,3) > 1,
    Fy = Fx(:,:,2);
    Fx = Fx(:,:,1);
end

fx = Fx - Fx([1 1:end-1],:);
fy = Fy - Fy(:,[1 1:end-1]);
fd = fx+fy;

end

