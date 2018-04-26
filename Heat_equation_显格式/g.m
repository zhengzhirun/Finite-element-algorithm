function z=g(p,t)
% 第一类边界条件
z = exp( 1/2 .* ( p(:,1) + p(:,2) ) - t );
end