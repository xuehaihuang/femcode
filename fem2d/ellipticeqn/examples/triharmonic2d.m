% This example is defined on the square domain (0, 1)\times (0,1)
clear all;
syms x y pi;
u=sin(pi*x)^3*sin(pi*y)^3;

gradu =  gradient(u, [x,y]);
laplaceu = laplacian(u, [x,y]);
laplace2u = laplacian(laplaceu, [x,y]);
laplace3u = laplacian(laplace2u, [x,y]);
f = - laplace3u;
f = simplify(f);

hessu = hessian(u,[x,y]);
uxxx = diff(u,x,3);
uyyy = diff(u,y,3);
uxxy = diff(hessu(1,1),y,1);
uxyy = diff(hessu(2,2),x,1);

ccode(f,'File','f.c');
ccode(u,'File','u.c');
ccode(gradu,'File','gradu.c');
ccode([hessu(1,1), hessu(2,2), hessu(1,2)],'File','hessu.c');
ccode([uxxx, uyyy, uxxy, uxyy],'File','grad3u.c');
% ccode(u)
% ccode(f)
% ccode(gradu)
% ccode([hessu(1,1), hessu(2,2), hessu(1,2)])
% ccode([uxxx, uyyy, uxxy, uxyy])

% f1=simplify(f11*2*mu+f12*lambda);
% f2=simplify(f21*2*mu+f22*lambda);