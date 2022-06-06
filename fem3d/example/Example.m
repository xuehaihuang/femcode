% Created by Xuehai Huang on Nov. 15, 2015
% Copyright 2015 WZU. All rights reserved.
% This example is defined on the cube domain (0, 1)\itmes(0,1)\itmes(0,1)
clear all;clc;
syms x y z mu lambda pii;
%mu=1;
%lambda=100;


X = [x y z];
w = sin(pii*x)^2*sin(pii*y)^2*sin(pii*z)^2;
W = w * [0;0;1];
u = curl(W,X);

gradu = jacobian(u,X);
gradu = simplify(gradu);

curlu = curl(u,X);
curlu = simplify(curlu);

gradcurlu = jacobian(curlu,X);
gradcurlu = simplify(gradcurlu);

f = curl(curlu,X);
f = simplify(f);

% write to files
ccode(f(1),'file','f1.c')
ccode(f(2),'file','f2.c')
ccode(f(3),'file','f3.c')
ccode(u(1),'file','u1.c')
ccode(u(2),'file','u2.c')
ccode(u(3),'file','u3.c')
ccode(gradu(1,1),'file','u1x.c')
ccode(gradu(1,2),'file','u1y.c')
ccode(gradu(1,3),'file','u1z.c')
ccode(gradu(2,1),'file','u2x.c')
ccode(gradu(2,2),'file','u2y.c')
ccode(gradu(2,3),'file','u2z.c')
ccode(gradu(3,1),'file','u3x.c')
ccode(gradu(3,2),'file','u3y.c')
ccode(gradu(3,3),'file','u3z.c')
ccode(gradu,'file','gradu.c')
ccode(curlu,'file','curlu.c')
ccode(gradcurlu,'file','gradcurlu.c')