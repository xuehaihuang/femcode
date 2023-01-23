% This example is defined on the square domain (-1, 1)\times (-1,1)
clear all;
syms x y mu lambda;
%mu=1;
%lambda=100;
w=sin(pi*x)^2*sin(pi*y)^2/2;

u1=diff(w,y,1);
u2=-diff(w,x,1);

u1x=diff(u1,x,1);
u1y=diff(u1,y,1);
u2x=diff(u2,x,1);
u2y=diff(u2,y,1);
eps11=u1x;
eps22=u2y;
eps12=(u1y+u2x)/2;
eps21=(u1y+u2x)/2;

traeps=eps11+eps22;

sigma11= 2*mu*eps11 + lambda*traeps;
sigma22= 2*mu*eps22 + lambda*traeps;
sigma12= 2*mu*eps12;
sigma21= 2*mu*eps21;

f11=diff(eps11,x,1)+diff(eps12,y,1);
f12=diff(traeps,x,1);
f21=diff(eps21,x,1)+diff(eps22,y,1);
f22=diff(traeps,y,1);

f11=-f11;
f12=-f12;
f21=-f21;
f22=-f22;

f1=simplify(f11*2*mu+f12*lambda);
f2=simplify(f21*2*mu+f22*lambda);