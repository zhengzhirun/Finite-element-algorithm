function show_energy(input_energy_file)
%% 文件读取操作
fid =  fopen(input_energy_file);
if fid < 0
    error('File %s not found',input_energy_file)
end

energy_cell = textscan(fid,'%f %f %f %f %f');
fclose(fid);

engmax = energy_cell{1};
engmin = energy_cell{2};
engavg = energy_cell{3};
engdev = energy_cell{4};
qavg = energy_cell{5};

N = size(engmax);
X = 1:N;

gcf = figure(1);
plot(X,engmax)
title('engmax');
saveas(gcf,'Engmax','jpg');

gcf = figure(2);
plot(X,engmin)
title('engmin');
saveas(gcf,'Engmin','jpg');

gcf = figure(3);
plot(X,engavg)
title('engavg');
saveas(gcf,'Engavg','jpg');

gcf = figure(4);
plot(X,engdev)
title('engdev')
saveas(gcf,'Engdev','jpg');

gcf = figure(5);
plot(X,qavg);
title('qavg')
saveas(gcf,'qavg','jpg');

