function z = exactu(p)   
% 解析解
z = exp(pi .* ( p(:,1) + p(:,2))) .* sin(pi.* p(:,1)) .* sin(pi.* p(:,2));
end
