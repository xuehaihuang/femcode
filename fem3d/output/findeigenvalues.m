clear all
load('A.dat')
%A=sparse(A(2:end,1),A(2:end,2),A(2:end,3),A(1,1),A(1,2),A(1,3));
A=sparse(A(2:end,1),A(2:end,2),A(2:end,3));
load('B.dat')
B=sparse(B(2:end,1),B(2:end,2),B(2:end,3));
num=min(4000, size(A, 1));
opts.tol = 1e-20;
E=eigs(A,B,num,'sm',opts);
E=sort(E);
F=E(E>0.005)
%F=eigs(A,B,100,'lm');
