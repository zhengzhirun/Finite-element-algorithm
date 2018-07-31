function u12 = crossProuct(v1,v2)
    N =size(v1,1);
    u12 = [zeros(N,1),zeros(N,1),v1(:,1).*v2(:,2)-v1(:,2).*v2(:,1)];
end
