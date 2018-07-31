function flag = genDensInfoH1(node,elem,ind)

%节点数目
N = size(node,1);
%单元数目
NT=size(elem,1);

%生成文件"dens_nodes.dat"
fid = fopen('dens_nodes.dat', 'w');

radi = 6371000;
x = 0.5*node(:,1)/radi;
y = 0.5*node(:,2)/radi;
xy2 = x.*x+y.*y;

sp_crd = zeros(N,3);
sp_crd(:,1) = radi*2*x./(1+xy2);
sp_crd(:,2) = radi*2*y./(1+xy2);
sp_crd(:,3) = radi*(1-xy2)./(1+xy2);
fprintf(fid, '%.12f  %.12f  %.12f\n', sp_crd'); 

fclose(fid);


%生成文件"dens_valus.dat"

%p2t(i,j)表示节点i在第j个单元里
p2t=sparse(elem,[1:NT,1:NT,1:NT],1,N,NT);

%生成单元尺寸(三边的平均长度)
ve1 = node(elem(:,3),:)-node(elem(:,2),:);
ve2 = node(elem(:,1),:)-node(elem(:,3),:);
ve3 = node(elem(:,2),:)-node(elem(:,1),:);
area = 0.5* sqrt(sum(crossProuct(ve3,-ve2).^2,2));
htri=sqrt(sum(ve1.^2,2))+sqrt(sum(ve2.^2,2))+sqrt(sum(ve3.^2,2));
htri=htri/3;

%计算rho
rho=(ind.^2)./(htri.^4);
rho=p2t*(rho.*area)./(p2t*area);
%h=figure;
for i=1:5
   rho=p2t'*rho/3;
   rho=p2t*(rho.*area)./(p2t*area);
end

fidd = fopen('dens_valus.dat', 'w');
fprintf(fidd, '%.12f\n', rho'); 
fclose(fidd);

flag = 1;

end
