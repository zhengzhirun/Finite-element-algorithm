clear all; close all; clc
%%

%参数设置
maxN = 2e3;
maxIt = 50;
theta = 0.5;
errL2 = zeros(maxIt,1);  %L2范
errH1 = zeros(maxIt,1);  %H1范

%生成初始网格
node = [0,0;0,1;1,0;1,1];  %点(全局编码)
elem = [3,4,1;2,1,4];   %边(局部编码)
elem = label(node,elem);  %给网格做上标签
gcf = figure(1);
showmesh(node,elem)
findelem(node,elem)
findnode(node,1:4)
title('初始网格')
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
pde.g_d = 0;

for k = 1:maxIt
    %step1 求解poisson问题
    u = poisson1(node,elem,pde);
    
    %画出网格
    gcf = figure(3);
    subplot(1,2,1) 
    showmesh(node,elem)
    %画出数值解
    subplot(1,2,2)
    showsolution(node,elem,u)
    pause(0.5)  %覆盖暂停0.5秒,观察细化过程
    
 %%   
    %step2 估计
    %eta = estimaterecovery(node,elem,u); %recovery type
    eta = estimateresidual(node,elem,u,pde); %residual type
    
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
    markedElem = mark(elem,eta,theta);   
%%   
    %step4 细化标记单元(二分法细化)
    [node,elem] = bisect(node,elem,markedElem);
end  
saveas(gcf,'数值结果','jpg')
%% 绘制收敛率曲线
N = N(1:k);
errH1 = errH1(1:k);
errL2 = errL2(1:k);
gcf = figure(4);
showrate2(N',errH1,3,'-*','||Du-Du_h||',N',errL2,3,'k-+','||u-u_h||');
title('收敛速率')
saveas(gcf,'收敛速率','jpg')   
    
    
    
    
    


