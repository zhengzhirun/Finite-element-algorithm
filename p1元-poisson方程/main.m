clear all; close all; clc
%%测试案例
pde.f = inline('-2 .* pi.^2 .* exp(pi .* (p(:,1) + p(:,2))) .* ( sin(pi.* p(:,1)).* cos(pi.* p(:,2)) + cos(pi.* p(:,1)) .* sin(pi.* p(:,2)) ) ','p');
pde.g_D = 0;

[node,elem] = squaremesh([0,1,0,1],0.03125);
figure(1)
gfc = showmesh(node,elem);
title('网格尺寸为dh=0.03125')
saveas(gcf,'网格剖分','jpg')

u = poisson(node,elem,pde);
N = size(node,1);

%作数值图
figure(2)
X = reshape(node(:,1),sqrt(N),sqrt(N));
Y = reshape(node(:,2),sqrt(N),sqrt(N));
U = reshape(u,sqrt(N),sqrt(N));
gcf = mesh(X,Y,U);
title('网格尺寸为dh=0.03125的数值解')
saveas(gcf,'数值解','jpg')

%做精确解的图
figure(3)
U1 = exp(pi .* ( X + Y)) .* sin(pi.*X) .* sin(pi.*Y);
gcf = mesh(X,Y,U1);
title('对应网格尺寸dh=0.03125的解析解')
saveas(gcf,'解析解','jpg')