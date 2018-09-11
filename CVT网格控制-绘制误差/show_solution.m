function show_solution(input_mesh_file,input_result_file)
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
title('实际计算的网格');
saveas(gfc,'实际计算的网格','jpg');

%% 文件读写操作(读入数据并且作图)
fid = fopen(input_result_file,'r');
if fid < 0
    error('File %s not found',input_result_file);
end
u_h_cell = textscan(fid,'%s');
fclose(fid);
u_h = u_h_cell{1};
u_h = str2double(u_h(:));
%% 作图
gcf = figure(2);
showsolution(node,elem,u_h);
axis([0 1 0 1 0 1]);
title('数值解结果');
saveas(gcf,'数值结算结果','jpg');

%%
gcf = figure(3);
u = u_exact(node);
showsolution(node,elem,u);
axis([0 1 0 1 0 1]);
title('解析解结果');
saveas(gcf,'解析解结果','jpg');