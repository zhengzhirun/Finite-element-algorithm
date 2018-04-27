clear all; close all; clc
%%
[node,elem] = squaremesh([0,1,0,1],1);
gcf = figure(1);
elem = label(node,elem);  %给网格做上标签
showmesh(node,elem)
findelem(node,elem)
findnode(node,1:4)
saveas(gcf,'初始网格','jpg')
errL2 = zeros(1,5);
errH1 = zeros(1,5);
N =zeros(1,5);

for k = 1:5
    %%
    [node,elem] = uniformrefine(node,elem);
    %% 参数
    pde.f = @f;
    pde.g_D = @g;
    pde.u0 = @u0;
    u = heat(node,elem,pde);  %求出结果
    %%
    errL2(k) = getL2error(node,elem,@exactu,u{1001});
    errH1(k) = getH1error(node,elem,@Du,u{1001});
    %%
    N(k) = size(node,1);
end
%%
gcf = figure(2);
showmesh(node,elem)
title('最终的计算网格')
saveas(gcf,'最终的计算网格','jpg')
%%
gcf = figure(3);
subplot(1,2,1)
showsolution(node,elem,u{1001})
title('时间0.1后的数值解')

subplot(1,2,2)
U = exp(1/2 .* (node(:,1)+node(:,2))-0.1);
showsolution(node,elem,U)
title('时间0.1后的解析解')

saveas(gcf,'数值解和解析解','jpg')
%%
gcf = figure(4);
N = N(1:k);
showrate2(N,errH1,3,'-*','||Du-Du_h||',N,errL2,3,'k-+','||u-u_h||');
title('时间t=0.1时的收敛曲线')
saveas(gcf,'收敛速率','jpg')

