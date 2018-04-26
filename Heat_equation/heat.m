function u = heat(node,elem,pde)
%%
%求解一个heat方程的有限元程序
dt = 0.001;
t = 0:0.001:1;
%
%  使用元胞数组来存储不同时间处的数值解U{0+dt}
u{1} = pde.u0(node);
for t_idx = 2:size(t,2)
    N = size(node,1);  NT = size(elem,1); 
    Ndof = N;

    %计算局部梯度算子
    [Dphi,area] = gradbasis(node,elem);

    %计算高斯积分点(阶为3)
    [lambda,weight] = quadpts(3);  
    phi = lambda;
    nQuad = size(lambda,1);
    
    %组装刚度矩阵A
    A = sparse(Ndof,Ndof);
    M = sparse(Ndof,Ndof);
    for i = 1:3
        for j = i:3
            % 计算a(u,v) 
            Aij = (Dphi(:,1,i).*Dphi(:,1,j) + Dphi(:,2,i).*Dphi(:,2,j)).*area;
        
            % 计算(u,v)
            Mij = zeros(NT,1);
            for p = 1:nQuad
                Mij = Mij + weight(p) * phi(p,i) * phi(p,j) .* area;
            end
            
            A = A + sparse(elem(:,i),elem(:,j),Aij,Ndof,Ndof);
            M = M + sparse(elem(:,i),elem(:,j),Mij,Ndof,Ndof);      
        end
    end
    clear K Aij Mij
    
    %% 组装右端项b
    b = zeros(Ndof,1);
    bt = zeros(NT,3);

    for j = 1:3
        for p = 1:nQuad
            pxy = lambda(p,1)*node(elem(:,1),:) ...
            + lambda(p,2)*node(elem(:,2),:) ...
            + lambda(p,3)*node(elem(:,3),:);
        bt(:,j) = bt(:,j) + weight(p) * pde.f(pxy,t(t_idx)) * phi(p,j);
        end
        bt(:,j) = bt(:,j) .* area;
    end
    b = accumarray(elem(:),bt(:),[Ndof 1]);
    clear pxy bt
    
    %变分形式出现的积分(类似Neumann边界条件处理)
	f = g_N(node,t(t_idx));
    if size(f,2)~=1
        f = f';   % transfer f to a column vector
    end
    if length(f)== N        % piecewise linear
       MF = accumarray([elem(:,1);elem(:,2);elem(:,3)],[area;area;area]/3,[N,1]);
           F = f.*MF;
    elseif length(f) == NT  % piecewise constant
           sf = f.*area;
           F = accumarray([elem(:,1);elem(:,2);elem(:,3)],[sf;sf;sf]/3,[N,1]);
    else
        error('The length of f should be either number of elements or number of nodes.')
    end
    b = b + F;
    

    %% 边界条件的处理(Dirichlet边界条件)(移到右端项,并对b进行修改)
    u{t_idx} = zeros(Ndof,1);
    fixedNode = []; %边界点 
    freeNode = [];  %内部点
    [fixedNode,bdEdge,isBdNode] = findboundary(elem);
    freeNode = find(~isBdNode);

    % AD(fixedNode,fixedNode)=I, AD(fixedNode,freeNode)=0, AD(freeNode,fixedNode)=0
    % MD(fixedNode,fixedNode)=I, MD(fixedNode,freeNode)=0, MD(freeNode,fixedNode)=0
    if ~isempty(fixedNode)
        bdidx = zeros(Ndof,1); 
        bdidx(fixedNode) = 1;
        Tbd = spdiags(bdidx,0,Ndof,Ndof);
        T = spdiags(1-bdidx,0,Ndof,Ndof);
        AD = T*A*T + Tbd;
        MD = T*A*T + Tbd;
    else
        AD = A;
        MD = M;
    end

    %% 修改右端项b
    u{t_idx}(fixedNode) = pde.g_D(node(fixedNode,:),t(t_idx));
    b = b * dt + MD * u{t_idx-1} - ( MD + AD * dt ) * u{t_idx};
    
    %
    %求解线性方程组
    AAD = MD + AD * dt;
    u{t_idx}(freeNode) = AAD(freeNode,freeNode)\b(freeNode);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%变分形式多出的一项,可以类似第二类边界条件的处理方式进行处理
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   function f = g_N(p,t)
%%
    f = zeros(size(p,1),1);
    x = p(:,1); y = p(:,2);
    uprime = [(exp(p(:,1)/2+p(:,2)/2 - t ) ./ 2)   (exp(p(:,1) ./ 2 + p(:,2) ./ 2 - t) ./ 2)];
    leftbd = (abs(x)<eps);  % n = (-1,0); 
    f(leftbd) = - uprime(leftbd,1);
    rightbd = (abs(x-1)<eps); % n = (1,0); 
    f(rightbd) = uprime(rightbd,1);
    topbd = (abs(y-1)<eps);   % n = (0,1)
    f(topbd) = uprime(topbd,2);
    bottombd = (abs(y)<eps);% n = (0,-1)
    f(bottombd) = - uprime(bottombd,2);
%%   
    end 
end