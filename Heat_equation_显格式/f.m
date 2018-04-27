function z=f(p,t)

% 右端项
z = -3/2 .* exp( 1/2 .* ( p(:,1) + p(:,2) ) -t );
end
