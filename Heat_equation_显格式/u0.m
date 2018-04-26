function z=u0(p)
%初值条件
z = exp( 1/2 .* ( p(:,1) + p(:,2) ) );
end