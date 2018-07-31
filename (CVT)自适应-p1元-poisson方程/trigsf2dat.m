function [ elem ] = trigsf2dat( datfile )
%trigSF2DAT read nodes from datfile

fid = fopen(datfile, 'r');
if fid < 0
  error('File %s not found', datfile);
end

%% read information of elems
tt = textscan(fid, '%d %d %d', 1);
NT = tt{1};
elems = textscan(fid, '%d %d %d %d', NT);
elem = {elems{2},elems{3},elems{4}};
elem = cell2mat(elem);


end

