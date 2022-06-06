clear all
clc
load A1.dat;load A2.dat;load A3.dat;load A4.dat;load A5.dat;load A6.dat;
load A7.dat;load A8.dat;load A9.dat;load A10.dat;load A11.dat;
load A12.dat;load dop.dat;
A1=sparse(A1(:,1),A1(:,2),A1(:,3));
A2=sparse(A2(:,1),A2(:,2),A2(:,3));
A3=sparse(A3(:,1),A3(:,2),A3(:,3));
A4=sparse(A4(:,1),A4(:,2),A4(:,3));
A5=sparse(A5(:,1),A5(:,2),A5(:,3));
A6=sparse(A6(:,1),A6(:,2),A6(:,3));
A7=sparse(A7(:,1),A7(:,2),A7(:,3));
A8=sparse(A8(:,1),A8(:,2),A8(:,3),size(A4,1),size(A1,1));
A9=sparse(A9(:,1),A9(:,2),A9(:,3));
A10=sparse(A10(:,1),A10(:,2),A10(:,3));
A11=sparse(A11(:,1),A11(:,2),A11(:,3),size(A4,1),size(A2,1));
A12=sparse(A12(:,1),A12(:,2),A12(:,3));
A1=full(A1);A2=full(A2);A3=full(A3);A4=full(A4);A5=full(A5);
A6=full(A6);A7=full(A7);A8=full(A8);A9=full(A9);A10=full(A10);
A11=full(A11);A12=full(A12);
D1=A1\A6';D2=A1\A7';D3=A1\A8';D4=A1\A9';
B1=A6*D1+A2;B2=A7*D2+A3;B3=A8*D3+A4;B4=A9*D4+A5;
B5=A7*D1+A10;B6=A8*D1+A11;B7=A9*D1+A12;
B8=A8*D2;B9=A9*D2;B10=A9*D3;
D5=B1\B5';D6=B1\B6';D7=B1\B7';
C1=B2-B5*D5;C2=B3-B6*D6;C3=B4-B7*D7;
C4=B8-B6*D5;C5=B9-B7*D5;C6=B10-B7*D6;

load BB1.dat;load BB2.dat;load BB3.dat;load BB4.dat;load BB5.dat;
load BB6.dat;load BB7.dat;load BB8.dat;load BB9.dat;load BB10.dat;
load CC1.dat;load CC2.dat;load CC3.dat;load CC4.dat;load CC5.dat;load CC6.dat;
load DD1.dat;load DD2.dat;load DD3.dat;load DD4.dat;load DD5.dat;load DD6.dat;load DD7.dat;
BB1=sparse(BB1(:,1),BB1(:,2),BB1(:,3));
BB2=sparse(BB2(:,1),BB2(:,2),BB2(:,3));
BB3=sparse(BB3(:,1),BB3(:,2),BB3(:,3));
BB4=sparse(BB4(:,1),BB4(:,2),BB4(:,3));
BB5=sparse(BB5(:,1),BB5(:,2),BB5(:,3));
BB6=sparse(BB6(:,1),BB6(:,2),BB6(:,3));
BB7=sparse(BB7(:,1),BB7(:,2),BB7(:,3));
BB8=sparse(BB8(:,1),BB8(:,2),BB8(:,3));
BB9=sparse(BB9(:,1),BB9(:,2),BB9(:,3));
BB10=sparse(BB10(:,1),BB10(:,2),BB10(:,3));
CC1=sparse(CC1(:,1),CC1(:,2),CC1(:,3));
CC2=sparse(CC2(:,1),CC2(:,2),CC2(:,3));
CC3=sparse(CC3(:,1),CC3(:,2),CC3(:,3));
CC4=sparse(CC4(:,1),CC4(:,2),CC4(:,3));
CC5=sparse(CC5(:,1),CC5(:,2),CC5(:,3));
CC6=sparse(CC6(:,1),CC6(:,2),CC6(:,3));
DD1=sparse(DD1(:,1),DD1(:,2),DD1(:,3));
DD2=sparse(DD2(:,1),DD2(:,2),DD2(:,3));
DD3=sparse(DD3(:,1),DD3(:,2),DD3(:,3));
DD4=sparse(DD4(:,1),DD4(:,2),DD4(:,3));
DD5=sparse(DD5(:,1),DD5(:,2),DD5(:,3));
DD6=sparse(DD6(:,1),DD6(:,2),DD6(:,3));
DD7=sparse(DD7(:,1),DD7(:,2),DD7(:,3));
%BB1=full(BB1);BB2=full(BB2);BB3=full(BB3);BB4=full(BB4);BB5=full(BB5);
%BB6=full(BB6);BB7=full(BB7);BB8=full(BB8);BB9=full(BB9);BB10=full(BB10);
CC1=full(CC1);CC2=full(CC2);CC3=full(CC3);CC4=full(CC4);CC5=full(CC5);CC6=full(CC6);
%DD1=full(DD1);DD2=full(DD2);DD3=full(DD3);DD4=full(DD4);DD5=full(DD5);DD6=full(DD6);DD7=full(DD7);
