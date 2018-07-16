%% 文件读写操作
%% 将网格信息写入txt文档中
clear all; close all; clc
[node,elem] = squaremesh([0,1,0,1],0.1);
figure(1)
gfc = showmesh(node,elem);

fid = fopen('Trimesh.txt','w');  % 以写的方式打开文件(文件不存在则自己创建)

%% 写入剖分节点
fprintf(fid,'%g\n',size(node,1));
for i = 1 : size(node,1)
    for j = 1 : 2
        fprintf(fid,'%g\t',node(i,j));
    end
    fprintf(fid,'\n');
end

%% 写入剖分单元(边)
fprintf(fid,'\n%g\n',size(elem,1));
for i = 1 : size(elem,1)
    for j = 1 : 3
        fprintf(fid,'%g\t',elem(i,j) - 1); % 和c++的编号统一(从0开始)
    end
    fprintf(fid,'\n');
end

%% 写入边界点
boundary_node = findboundary(elem);
boundary_node = boundary_node - 1;  % 和c++的编号统一(从0开始)
fprintf(fid,'\n%g\n',size(boundary_node,1));
for i = 1 : size(boundary_node,1)
    fprintf(fid,'%g\n',boundary_node(i));
end

fclose(fid);


