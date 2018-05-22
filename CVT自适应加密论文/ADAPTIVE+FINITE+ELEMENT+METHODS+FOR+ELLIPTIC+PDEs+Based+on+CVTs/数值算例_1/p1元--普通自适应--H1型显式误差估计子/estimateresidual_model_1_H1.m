function eta = estimateresidual_model_1_H1(node,elem,u,pde)
%%模型问题1的H_1型的误差估计子的计算


%% 计算有限元函数u的梯度
[Du,area] = gradu(node,elem,u);
center = (node(elem(:,1),:) + node(elem(:,2),:) + node(elem(:,3),:))/3;
K = a(center);              
Du = [K.*Du(:,1),K.*Du(:,2)];   

%% 计算跳量
% data structure
T = auxstructure(elem);
neighbor = T.neighbor; 
clear T
% edge vector
ve(:,:,1) = node(elem(:,3),:)-node(elem(:,2),:);
ve(:,:,2) = node(elem(:,1),:)-node(elem(:,3),:);
ve(:,:,3) = node(elem(:,2),:)-node(elem(:,1),:);
% scaled normal vector: a right 90 degree rotation of edge vector
ne(:,1,1)= ve(:,2,1); ne(:,2,1)= -ve(:,1,1);
ne(:,1,2)= ve(:,2,2); ne(:,2,2)= -ve(:,1,2);
ne(:,1,3)= ve(:,2,3); ne(:,2,3)= -ve(:,1,3);
clear ve;
% for boundary edges e, neighbor(t,e) = t. So the difference is zero.     
edgeJump = dot((Du-Du(neighbor(:,1),:)),ne(:,:,1),2).^2 ...
         + dot((Du-Du(neighbor(:,2),:)),ne(:,:,2),2).^2 ...
         + dot((Du-Du(neighbor(:,3),:)),ne(:,:,3),2).^2;

%% Elementwise residual
elemResidual = zeros(size(elem,1),1);


[lambda,weight] = quadpts(3);
nQuad = size(lambda,1);
for p = 1:nQuad
    % quadrature points in the x-y coordinate
    pxy = lambda(p,1)*node(elem(:,1),:) ...
        + lambda(p,2)*node(elem(:,2),:) ...
        + lambda(p,3)*node(elem(:,3),:);
    fp = pde.f(pxy);
    
    bp = b(pxy);   %b是一个函数

    uhp = u(elem(:,1))*lambda(p,1) + ...
			u(elem(:,2))*lambda(p,2) + ...
			u(elem(:,3))*lambda(p,3);

    zp = fp - bp .* uhp  - 10 .* sin(pxy(:,2)) .* Du(:,2) ;
    elemResidual = elemResidual + weight(p) * zp.^2 ;
end
elemResidual = elemResidual.*(area.^2);

%% Residual type error estimator
eta = (elemResidual + edgeJump).^(1/2);