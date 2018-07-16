%% 根据真解求函数的右端项f和pde方程的边界条件g
clear all; close all; clc
syms x y
u = 1.0 / ((x-0.5)^2+(y-0.5)^2+0.01) - 1.0 / ((x+0.5)^2+(y+0.5)^2+0.01);
a = 10.0 * cos(y);
b = x^2 + y^2;

dadx = diff(a,x,1);
dady = diff(a,y,1);
dudx = diff(u,x,1);
dudy = diff(u,y,1);

dudx2 = diff(u,x,2);
dudy2 = diff(u,y,2);

%右端项
f = -(dadx * dudx + a * dudx2 + dady * dudy + a * dudy2) + b * u;
f = simplify(f)

%边值条件
g = u;
g = simplify(g) 