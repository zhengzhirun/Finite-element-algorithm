clear all; clc; close all;
%% 参数设置
maxN = 2e3;
maxIt = 50;
theta = 0.5;
errL2 = zeros(maxIt,1);  %L2范
errH1 = zeros(maxIt,1);  %H1范
%% 初始化网格
system('./mesh_init -f "square.poly" -q 1 -a 0.2 ');
% 得到网格信息 
[ node ] = nodesf2dat( 'nodes.dat' ); 
[ elem ] = trigsf2dat( 'trigs.dat' );
elem = double(elem);
%% 显示初始网格(存储初始网格)
showmesh(node,elem)
saveas(gcf,'初始网格','jpg')
%%
%问题边界条件以及右端项
pde.f = @f;
pde.g_d = 0;

for k = 1:maxIt
    %step1 求解poisson问题
    u = poisson1(node,elem,pde);
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
    %% Step 3: 优化和细化网格
    genDensInfoH1(node,elem,eta);
    system('./mesh_refi -b 1 -p 0.2');
    system('./mesh_opti -b 1');
    [ node ] = nodesf2dat( 'nodes.dat' );
	[ elem ] = trigsf2dat( 'trigs.dat' );
	elem = double(elem);
end
%% 画出最终网格
gcf = figure(3);
subplot(1,2,1) 
showmesh(node,elem)
%画出数值解
subplot(1,2,2)
showsolution(node,elem,u)
saveas(gcf,'数值结果','jpg')

%% 绘制收敛率曲线
N = N(1:k);
errH1 = errH1(1:k);
errL2 = errL2(1:k);
gcf = figure(4);
showrate2(N',errH1,3,'-*','||Du-Du_h||',N',errL2,3,'k-+','||u-u_h||');
title('收敛速率')
saveas(gcf,'收敛速率','jpg')   
    




