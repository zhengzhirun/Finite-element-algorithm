%% 网格信息
clear all; close all; clc

[node,elem] = squaremesh([0,1,0,1],1/64);
figure(1)
gfc = showmesh(node,elem);
title('实际计算的网格');
saveas(gfc,'实际计算的网格','jpg');

%% 文件读写操作(读入数据并且作图)
fid = fopen('Results.txt');
u_h_cell = textscan(fid,'%s');
fclose(fid);
u_h = u_h_cell{1};
u_h = str2double(u_h(:));

%% 作图
gcf = figure(2);
showsolution(node,elem,u_h);
title('数值解结果');
axis([0 1 0 1 0 60]);
saveas(gcf,'数值结算结果','jpg');
%%
gcf = figure(3);
u = u_exact(node);
showsolution(node,elem,u);
title('解析解结果');
axis([0 1 0 1 0 60]);
saveas(gcf,'解析解结果','jpg');