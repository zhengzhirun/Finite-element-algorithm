function value = a(p)
x = p(:,1); y = p(:,2);
value = 10.0 .* cos(y) + 0.*x;
end