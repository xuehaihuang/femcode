% This example is defined on the square domain (0, 1)\times (0,1)
clear all;clc;
syms x y z mu lambda pi;
%mu=1;
%lambda=100;


X = [x y z];
w = sin(pi*x)^2*sin(pi*y)^2*sin(pi*z)^2 / pi;
W = w * [0;0;1];
u = curl(W,X);

gradu = jacobian(u,X);
gradu = simplify(gradu);

laplaceu = [laplacian(u(1), X);laplacian(u(2), X);laplacian(u(3), X)];

f = -laplaceu;
f = simplify(f);

% write to files
ccode(f,'file','f.c')
ccode(u,'file','u.c')
ccode(gradu,'file','gradu.c')
% ccode(gradu(1,1),'file','u1x.c')
% ccode(gradu(1,2),'file','u1y.c')
% ccode(gradu(1,3),'file','u1z.c')
% ccode(gradu(2,1),'file','u2x.c')
% ccode(gradu(2,2),'file','u2y.c')
% ccode(gradu(2,3),'file','u2z.c')
% ccode(gradu(3,1),'file','u3x.c')
% ccode(gradu(3,2),'file','u3y.c')
% ccode(gradu(3,3),'file','u3z.c')
% ccode(gradu,'file','gradu.c')
% ccode(curlu,'file','curlu.c')
% ccode(gradcurlu,'file','gradcurlu.c')