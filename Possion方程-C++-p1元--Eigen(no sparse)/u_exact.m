function value = u_exact(node)
    value = exp( (pi .* ( node(:,1) + node(:,2) ))) .* sin(pi .* node(:,1)) .* sin(pi .* node(:,2));
end