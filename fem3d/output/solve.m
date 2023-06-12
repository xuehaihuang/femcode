clear all
clc
load A.dat;
load b.dat;
i=A(:,1);
j=A(:,2);
s=A(:,3);
clear A;
A=sparse(i,j,s);
clear i j s;
u=A\b;
dlmwrite('uh.dat',u,'\n')

