clear;

h=1:8;
h=(1/2).^h;

errorsigma0=[3.7996E-02	3.2925E-02	2.0099E-02	1.0639E-02	5.4273E-03	2.7371E-03	1.3741E-03	6.8838E-04];
erroru1=[4.3529E-03	2.6350E-03	9.6374E-04	2.7287E-04	7.1159E-05	1.8078E-05	4.5499E-06	1.1372E-06];
erroruE=[5.5396E-02	3.6198E-02	1.9442E-02	9.7570E-03	4.8673E-03	2.4302E-03	1.2143E-03	6.0694E-04];


h=-log(h);
errorsigma0=log(errorsigma0);
erroru1=log(erroru1);
erroruE=log(erroruE);

plot(h,errorsigma0,'o-',h,erroru1,'<-',h,erroruE,'*-','LineWidth',2);
leg=legend('$\ln\|\mathbf{\sigma}-\mathbf{\sigma}_h\|_0$','$\ln|u-u_h|_{1}$','$\ln|\!|\!|u-u_h|\!|\!|$');
%leg=legend('$\|\mathbf{\sigma}-\mathbf{\sigma}_h\|_0$','$|u-u_h|_{1}$','$|\!|\!|u-u_h|\!|\!|$','Location','southwest');
set(leg,'Interpreter','latex');
set(xlabel('$\ln(1/h)$'),'Interpreter','latex');
%set(ylabel('$\ln$Errors'),'Interpreter','latex');

%set (gcf,'Position',[200, 100, 720, 540])
%set (gcf,'Position',[200, 100, 800, 600], 'color','w')
% reset: set (gcf,'Position',[232, 246, 560, 420], 'color','w')


line([1 1 2 1]+3,[-4 -5 -5 -4]-2.5,'LineWidth',2,'Color',[0 0 0]);

line([1 1 2 1]+2.5,[-7 -9 -9 -7]-3.5,'LineWidth',2,'Color',[0 0 0]);
%line([3.5 2.5 3.5 3.5],[-3 -3 -4 -3],'LineWidth',2,'Color',[0 0 0]);

%line([2.5 2.5 3.5 2.5],[-1 -2 -2 -1],'LineWidth',2,'Color',[0 0 0]);

%line([6 6 8 6],[0 -1 -1 0],'LineWidth',2,'Color',[0 0 0]);



%plot(h,errorEk3l0,'p-',h,errorEk3l1,'<-',h,errorEk3l2,'o-',h,errorEk3l3,'s-','LineWidth',2);
%line([1.5 1.5 2.5 1.5],[-6 -8 -8 -6],'LineWidth',2,'Color',[0 0 0]);
%legend('l=0', 'l=1', 'l=2', 'l=3');

%subplot(2,1,1)
%plot(h,error1k2l1,'p-',h,error1k2l0,'<-',h,error1k1l0,'o-','LineWidth',2);
%legend('k=2,l=1', 'k=2,l=0', 'k=1,l=0', 'Location', 'SouthWest');
%grid on;
%subplot(2,1,2)
%plot(h,errorEk2l1,'p-',h,errorEk2l0,'<-',h,errorEk1l0,'o-','LineWidth',2);
%legend('k=2,l=1', 'k=2,l=0', 'k=1,l=0', 'Location', 'SouthWest');

grid on;

