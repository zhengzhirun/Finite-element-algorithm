function show_mesh(input_mesh_file,i)
%% 文件读入操作
fid = fopen(input_mesh_file,'r');
if fid < 0
    error('File %s not found',input_mesh_file);
end
np = textscan(fid,'%d',1);
N = np{1};
points = textscan(fid,'%f %f',N);
points = {points{1},points{2}};
node = cell2mat(points);
%%
mp = textscan(fid,'%d,%d,%d,%d',N+1);
M = mp{1};
elements = textscan(fid,'%d %d %d',M);
elements = {elements{1},elements{2},elements{3}};
elem = cell2mat(elements);
elem = elem + 1;
fclose(fid);
%%
figure(1)
gfc = showmesh(node,elem);
index = int2str(i);
title(index);
saveas(gfc,index,'jpg');


