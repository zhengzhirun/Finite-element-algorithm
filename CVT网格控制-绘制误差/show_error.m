function show_error(input_error_file)
%% 文件读取操作
fid = fopen(input_error_file);
if fid < 0
    error('File %s not found',input_mesh_file);
end

Error_cell = textscan(fid,'%f %f %f %f');
fclose(fid);
L2Error = Error_cell{2};
H1Error = Error_cell{3};
index = Error_cell{1};

gcf = figure(1);
plot(index,L2Error)
title('L2Error')
saveas(gcf,'L2Error','jpg');

gcf = figure(2);
plot(index,H1Error)
title('H1Error')
saveas(gcf,'H1Error','jpg');
