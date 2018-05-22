clear all; close all; clc
%%测试案例
pde.f = @f;
pde.g_D = @g;

[node,elem] = squaremesh([-1,1,-1,1],1);
for i = 1:1
    [node,elem] = uniformbisect(node,elem);
end

figure(1)
gfc = showmesh(node,elem);
title('初始网格')
saveas(gcf,'初始计算网格','jpg')
%%
errL2 = zeros(1,7);
errH1 = zeros(1,7);
N = zeros(1,7);
for i = 1:7
    %对模型问题1进行计算
    [node,elem] = uniformbisect(node,elem);
    u = model_1(node,elem,pde);
    errL2(i) = getL2error(node,elem,@exactu,u);
    errH1(i) = getH1error(node,elem,@Du,u);
    N(i) = size(node,1);
end

%作出数值解的图
gcf = figure(2);
subplot(1,2,1)
showsolution(node,elem,u)
title('数值解')

%作出精确解的图
subplot(1,2,2)
U = exactu(node);
showsolution(node,elem,U)
title('解析解')
saveas(gcf,'数值解和解析解','jpg')
%%
%画出收敛曲线图
gcf = figure(3);
N = N(1:i);
showrate2(N,errH1,3,'-*','||Du-Du_h||',N,errL2,3,'k-+','||u-u_h||');
saveas(gcf,'收敛速率','jpg')
