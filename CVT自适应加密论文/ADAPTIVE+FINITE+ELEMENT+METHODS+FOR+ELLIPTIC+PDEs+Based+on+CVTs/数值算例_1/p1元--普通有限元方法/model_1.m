function u = model_1(node,elem,pde)
%求解模型问题1的有限元程序
N = size(node,1);  NT = size(elem,1); 
Ndof = N;

%计算梯度算子和三角单元的面积
[Dphi,area] = gradbasis(node,elem);

A = sparse(N,N);  %刚度矩阵A使用稀疏矩阵存储方式
B = sparse(N,N);
D = sparse(N,N);
%组装刚度矩阵A
for i = 1:3
    for j = 1:3
        %% 第一项的计算
        
        % compute quadrature points
        [lambda,weight] = quadpts(3);
        nQuad = size(lambda,1);
        % compute element-wise coefficient d
        dt = zeros(NT,size(a([0 0]),2));
        for p = 1:nQuad
        % quadrature points in the x-y coordinate
            pxy = lambda(p,1)*node(elem(:,1),:) ...
                + lambda(p,2)*node(elem(:,2),:) ...
                + lambda(p,3)*node(elem(:,3),:);
            dt = dt + weight(p)*a(pxy);
            if size(dt,2) == 1     % scalar function
                Dij = dt.*(Dphi(:,1,i).*Dphi(:,1,j) + Dphi(:,2,i).*Dphi(:,2,j));
            end
        end
        Dij = Dij .* area;
        D = D + sparse(elem(:,i),elem(:,j),Dij,N,N);
        
        %%第二项
        % compute quadrature points
        [lambda,weight] = quadpts(3);
        phi = lambda;
        nQuad = size(lambda,1);
        % compute element-wise coefficient B
        Bij = zeros(NT,1);
        for p = 1:nQuad
            % quadrature points in the x-y coordinate
            pxy = lambda(p,1)*node(elem(:,1),:) ...
                    + lambda(p,2)*node(elem(:,2),:) ...
                    + lambda(p,3)*node(elem(:,3),:);
            Bij = Bij + weight(p)*phi(p,i)*phi(p,j)*b(pxy);
        end
        Bij = Bij .* area;
        B = B + sparse(elem(:,i),elem(:,j),Bij,N,N);
          
    end  %end for j
end %end for i


% 组装右端项F
F = zeros(Ndof,1);
[lambda,weight] = quadpts(3);  %参见quadpts.html
phi = lambda;
nQuad = size(lambda,1);
bt = zeros(NT,3);
for p = 1:nQuad
    % quadrature points in the x-y coordinate
    pxy = lambda(p,1)*node(elem(:,1),:) ...
        + lambda(p,2)*node(elem(:,2),:) ...
        + lambda(p,3)*node(elem(:,3),:);
    fp = pde.f(pxy);
    for i = 1:3
        bt(:,i) = bt(:,i) + weight(p)*phi(p,i)*fp;
    end
end
    bt = bt.*repmat(area,1,3);
    F = accumarray(elem(:),bt(:),[Ndof 1]);
clear pxy bt

% 边界条件的处理(移到右端项,并对b进行修改)
u = zeros(Ndof,1);
fixedNode = []; %边界点 
freeNode = [];  %内部点
[fixedNode,bdEdge,isBdNode] = findboundary(elem);
freeNode = find(~isBdNode);

% AD(fixedNode,fixedNode)=I, AD(fixedNode,freeNode)=0, AD(freeNode,fixedNode)=0
A = D + B;
if ~isempty(fixedNode)
    bdidx = zeros(Ndof,1); 
    bdidx(fixedNode) = 1;
    Tbd = spdiags(bdidx,0,Ndof,Ndof);
    T = spdiags(1-bdidx,0,Ndof,Ndof);
    AD = T*A*T + Tbd;
else
    AD = A;
end

% 修改右端项F
u(fixedNode) = pde.g_D(node(fixedNode,:));  %第一类边界条件,边界处的值已经给定
F = F - A*u;

%求解线性方程组
u(freeNode) = AD(freeNode,freeNode)\F(freeNode);
end