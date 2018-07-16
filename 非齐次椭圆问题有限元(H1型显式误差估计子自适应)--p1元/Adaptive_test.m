clear all; close all; clc
%%

%参数设置
maxN = 2e4;
maxIt = 50;
theta = 0.5;
errL2 = zeros(maxIt,1);  %L2范
errH1 = zeros(maxIt,1);  %H1范

%生成初始网格
[node,elem] = squaremesh([-1,1,-1,1],1);
showmesh(node,elem)
gcf = figure(1);
showmesh(node,elem)
title('初始网格');
saveas(gcf,'初始网格','jpg')
%%
%通过对网格一致二分加密得到细网格(两种加密形式,任选一种)
for k =1:2
    [node,elem] = uniformbisect(node,elem);
    %[node,elem] = uniformrefine(node,elem);   %加密网格
end
gcf = figure(2);
showmesh(node,elem)
title('初始网格一致二分加密两次')
saveas(gcf,'初始网格一致二分加密两次','jpg')
%%
%问题边界条件以及右端项
pde.f = @f;
pde.g_D = @g;

for k = 1:maxIt
    %step1 求解poisson问题
    u = model_1(node,elem,pde);
%%
    %step2 估计
    %eta = estimaterecovery(node,elem,u); %recovery type
    eta = estimateresidual_model_1_H1(node,elem,u,pde); %residual type
    
    %记录误差和点的数量
    errL2(k) = getL2error(node,elem,@exactu,u);
    errH1(k) = getH1error(node,elem,@Du,u);
    N(k) = size(node,1);
    if N(k) > maxN
       break;
    end
%%    
    %step3 标记细化单元
    %默认'L2'加密,对残差前50%加密(从到小到大排序后的前50%)
    if k >= maxIt
        break
    end
        markedElem = mark(elem,eta,theta);   
%%   
    %step4 细化标记单元(二分法细化)
    [node,elem] = bisect(node,elem,markedElem);
    figure(6)
    showmesh(node,elem)
    pause(0.5)
end
%%
%% 画出最后计算结果的数值解和解析解
    gcf = figure(3);
    subplot(1,2,1)
    showsolution(node,elem,u)
    title('数值解')

    subplot(1,2,2)
    U = exactu(node);
    showsolution(node,elem,U)
    title('解析解')
    pause(0.5)
    saveas(gcf,'数值结果','jpg')
%% 画出最后的计算网格
gcf = figure(4);
showmesh(node,elem)
title('最终的计算网格')
saveas(gcf,'最终的计算网格','jpg')


%% 绘制收敛率曲线
N = N(1:k);
errH1 = errH1(1:k);
errL2 = errL2(1:k);
gcf = figure(5);
showrate2(N',errH1,2,'-*','||Du-Du_h||',N',errL2,2,'k-+','||u-u_h||');
title('收敛速率')
saveas(gcf,'收敛速率','jpg')   