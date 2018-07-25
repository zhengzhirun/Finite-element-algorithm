%% 文件读写操作
%% 将网格信息写入txt文档中
clear all; close all; clc
[node,elem] = squaremesh([0,1,0,1],0.1);

fid = fopen('Trimesh.txt','w');  % 以写的方式打开文件(文件不存在则自己创建)
%% 写入剖分节点
fprintf(fid,'%d\n',size(node,1));
fprintf(fid,'%g\t%g\n',node');

%% 写入剖分单元(边)
fprintf(fid,'\n%g\n',size(elem,1));
fprintf(fid,'%g\t%g\t%g\n',(elem - 1)');

%% 写入边界点
boundary_node = findboundary(elem);
boundary_node = boundary_node - 1;  % 和c++的编号统一(从0开始)
fprintf(fid,'\n%g\n',size(boundary_node,1));
%%
fprintf(fid,'%g\n',boundary_node');
fclose(fid);