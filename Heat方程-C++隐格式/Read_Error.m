%% 网格信息
clear all; close all; clc
filename{1} = 'Error_0.txt';
filename{2} = 'Error_1.txt';
filename{3} = 'Error_2.txt';
filename{4} = 'Error_3.txt';
filename{5} = 'Error_4.txt';
for i = 1:5
    fid = fopen(filename{i});
    Error_cell = textscan(fid,'%s %s %s');
    fclose(fid);
    L2Error = Error_cell{2}(2:end);
    L2Error = str2double(L2Error(:));
    L2ErrorNew(i) = L2Error(end);
    H1Error = Error_cell{3}(2:end);
    H1Error = str2double(H1Error(:));
    H1ErrorNew(i) = H1Error(end);
end

Nodes = [9 25 81 289 1089];
gcf = figure(1);
N = Nodes(1:5);
showrate2(N,H1ErrorNew,2,'-*','||Du-Du_h||',N,L2ErrorNew,2,'k-+','||u-u_h||');
title('时间t=1时的收敛曲线')
saveas(gcf,'收敛速率','jpg')