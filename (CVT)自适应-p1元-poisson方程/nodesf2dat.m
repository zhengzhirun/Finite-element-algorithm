function [ node ] = nodesf2dat( datfile )
%NODESF2DAT read nodes from datfile

fid = fopen(datfile, 'r');
if fid < 0
  error('File %s not found', datfile);
end

%% read information of points
np = textscan(fid, '%d %d %d %d', 1);
N = np{1};
points = textscan(fid, '%d %f %f %d', N);
points = {points{2},points{3}};
node = cell2mat(points);

end

